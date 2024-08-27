# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Quaternion provides a class for manipulating quaternion objects.  This class provides:

   - a convenient constructor to convert to/from Euler Angles (RA,Dec,Roll)
       to/from quaternions
   - class methods to multiply and divide quaternions

:Copyright: Smithsonian Astrophysical Observatory (2010)
:Author: Jean Connelly (jconnelly@cfa.harvard.edu)
"""
## Copyright (c) 2010, Smithsonian Astrophysical Observatory
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.
##     * Redistributions in binary form must reproduce the above copyright
##       notice, this list of conditions and the following disclaimer in the
##       documentation and/or other materials provided with the distribution.
##     * Neither the name of the Smithsonian Astrophysical Observatory nor the
##       names of its contributors may be used to endorse or promote products
##       derived from this software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
## ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
## WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
## DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
## (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
## LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
## ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
## SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import numpy as np
import warnings


class Quat(object):
    """
     Quaternion class

     Example usage::

      >>> from Quaternion import Quat
      >>> quat = Quat((12,45,45))
      >>> quat.ra, quat.dec, quat.roll
      (12, 45, 45)
      >>> quat.q
      array([ 0.38857298, -0.3146602 ,  0.23486498,  0.8335697 ])
      >>> q2 = Quat(quat.q)
      >>> q2.ra
      12.0

     Multiplication and division operators are overloaded for the class to
     perform appropriate quaternion multiplication and division.

     Quaternion composition as a multiplication q = q1 * q2 is equivalent to
     applying the q2 transform followed by the q1 transform.  Another way to
     express this is::

      q = Quat(numpy.dot(q1.transform, q2.transform))

     Example usage::

      >>> q1 = Quat((20, 30, 0))
      >>> q2 = Quat((0, 0, 40))
      >>> (q1 * q2).equatorial
      array([20., 30., 40.])

    This example first rolls about X by 40 degrees, then rotates that rolled
    frame to RA=20 and Dec=30.  Doing the composition in the other order does
    a roll about (the original) X-axis of the (RA, Dec) = (20, 30) frame,
    yielding a non-intuitive though correct result::

      >>> (q2 * q1).equatorial
      array([ 353.37684725,   34.98868888,   47.499696  ])


    :param attitude: initialization attitude for quat

    ``attitude`` may be:
      * another Quat
      * a 4 element array (expects x,y,z,w quat form)
      * a 3 element array (expects ra,dec,roll in degrees)
      * a 3x3 transform/rotation matrix
     N x those types :
      * an N x 4 element array ( N by x,y,z,w quat form)
      * an N x 3 element array ( N by ra, dec, roll in degrees,
        if N == 3, the optional 'intype = 'equatorial' may be used to differentiate,
        this from a transform matrix)
      * an N x 3x3

    :param intype: optional type to describe input attitude

    ``intype`` may be:
      * transform
      * equatorial
      * quaternion

    """

    def __init__(self, attitude, intype=None):
        warnings.warn(
            "mica.quaternion is deprecated as of Quaternion 3.5.0", FutureWarning
        )
        self._q = None
        self._equatorial = None
        self._ra0 = None
        self._roll0 = None
        self._T = None
        # checks to see if we've been passed a Quat
        if isinstance(attitude, Quat):
            self._set_q(attitude.q)
        else:
            # make it an array and check to see if it is a supported shape
            attitude = np.array(attitude)
            if (
                attitude.shape == (3, 3) and (intype is None or intype == "transform")
            ) or (
                attitude.ndim == 3
                and attitude.shape[-1] == 3
                and attitude.shape[-2] == 3
            ):
                self._set_transform(attitude)
            elif (
                intype == "quaternion"
                or attitude.shape == (4,)
                or (attitude.ndim == 2 and attitude.shape[-1] == 4)
            ):
                self._set_q(attitude)
            elif (
                intype == "equatorial"
                or attitude.shape == (3,)
                or (attitude.ndim == 2 and attitude.shape[-1] == 3)
            ):
                self._set_equatorial(attitude)
            else:
                raise TypeError(
                    "attitude is not one of possible types"
                    " (3 or 4 elements, Quat, or 3x3 matrix or N x (those types))"
                )

    def _set_q(self, q):
        """
        Set the value of the 4 element quaternion vector
        May be 4 element list or array or N x 4 element array.

        :param q: list or array of normalized quaternion elements

        """
        q = np.array(q)
        if q.ndim == 1:
            q = q[np.newaxis]
        if np.any((np.sum(q * q, axis=-1)[:, np.newaxis] - 1.0) > 1e-6):
            raise ValueError(
                "Quaternions must be normalized so sum(q**2) == 1;"
                " use Quaternion.normalize"
            )
        self._q = q
        flip_q = q[:, 3] < 0
        self._q[flip_q] = -1 * q[flip_q]
        # Erase internal values of other representations
        self._equatorial = None
        self._ra0 = None
        self._roll0 = None
        self._T = None

    def _get_q(self):
        """
        Retrieve 4-vector of quaternion elements in [x, y, z, w] form
        or N x 4-vector if N > 1.

        :rtype: numpy array

        """
        if self._q is None:
            # Figure out q from available values
            if self._equatorial is not None:
                self._q = self._equatorial2quat()
            elif self._T is not None:
                self._q = self._transform2quat()
        if self._q.shape[0] == 1:
            return self._q[0, :]
        return self._q

    # use property to make this get/set automatic
    q = property(_get_q, _set_q)

    def _set_equatorial(self, equatorial):
        """
        Set the value of the 3 element equatorial coordinate list [RA,Dec,Roll]
        expects values in degrees
        bounds are not checked

        :param equatorial: list or array [ RA, Dec, Roll] in degrees

        """
        att = np.array(equatorial)
        if att.ndim == 1:
            att = att[np.newaxis]
        ra, dec, roll = att[:, 0], att[:, 1], att[:, 2]
        self._ra0 = ra
        self._roll0 = roll
        self._ra0[ra > 180] = ra - 360
        self._roll0[roll > 180] = roll - 360
        self._equatorial = att

    def _get_equatorial(self):
        """Retrieve [RA, Dec, Roll]

        :rtype: numpy array
        """
        if self._equatorial is None:
            if self._q is not None:
                self._equatorial = self._quat2equatorial()
            elif self._T is not None:
                self._q = self._transform2quat()
                self._equatorial = self._quat2equatorial()
        if self._equatorial.shape[0] == 1:
            return self._equatorial[0, :]
        return self._equatorial

    equatorial = property(_get_equatorial, _set_equatorial)

    def _get_ra(self):
        """Retrieve RA term from equatorial system in degrees"""
        if self._equatorial is None:
            self.equatorial
        if self._equatorial.shape[0] == 1:
            return self._equatorial[:, 0][0]
        else:
            return self._equatorial[:, 0]

    def _get_dec(self):
        """Retrieve Dec term from equatorial system in degrees"""
        if self._equatorial is None:
            self.equatorial
        if self._equatorial.shape[0] == 1:
            return self._equatorial[:, 1][0]
        else:
            return self._equatorial[:, 1]

    def _get_roll(self):
        """Retrieve Roll term from equatorial system in degrees"""
        if self._equatorial is None:
            self.equatorial
        if self._equatorial.shape[0] == 1:
            return self._equatorial[:, 2][0]
        else:
            return self._equatorial[:, 2]

    ra = property(_get_ra)
    dec = property(_get_dec)
    roll = property(_get_roll)

    def _set_transform(self, T):
        """
        Set the value of the 3x3 rotation/transform matrix

        :param T: 3x3 array/numpy array
        """
        transform = np.array(T)
        if transform.ndim == 2:
            transform = transform[np.newaxis]
        self._T = transform

    def _get_transform(self):
        """
        Retrieve the value of the 3x3 rotation/transform matrix

        :returns: 3x3 rotation/transform matrix
        :rtype: numpy array

        """
        if self._T is None:
            if self._q is not None:
                self._T = self._quat2transform()
            elif self._equatorial is not None:
                self._T = self._equatorial2transform()
        if self._T.shape[0] == 1:
            return self._T[0, :]
        return self._T

    transform = property(_get_transform, _set_transform)

    def _quat2equatorial(self):
        """
        Determine Right Ascension, Declination, and Roll for the
        object quaternion

        :returns: N x (RA, Dec, Roll)
        :rtype: numpy array [ra,dec,roll]
        """

        q = self.q
        if q.ndim == 1:
            q = q[np.newaxis]
        q2 = q**2

        ## calculate direction cosine matrix elements from $quaternions
        xa = q2[:, 0] - q2[:, 1] - q2[:, 2] + q2[:, 3]
        xb = 2 * (q[:, 0] * q[:, 1] + q[:, 2] * q[:, 3])
        xn = 2 * (q[:, 0] * q[:, 2] - q[:, 1] * q[:, 3])
        yn = 2 * (q[:, 1] * q[:, 2] + q[:, 0] * q[:, 3])
        zn = q2[:, 3] + q2[:, 2] - q2[:, 0] - q2[:, 1]

        ##; calculate RA, Dec, Roll from cosine matrix elements
        ra = np.degrees(np.arctan2(xb, xa))
        dec = np.degrees(np.arctan2(xn, np.sqrt(1 - xn**2)))
        roll = np.degrees(np.arctan2(yn, zn))

        ra[ra < 0] = ra[ra < 0] + 360
        roll[roll < 0] = roll[roll < 0] + 360
        return np.array([ra, dec, roll]).transpose()

    #  _quat2transform is largely from Enthought's quaternion.rotmat, though this math is
    #  probably from Hamilton.
    #  License included for completeness
    #
    # This software is OSI Certified Open Source Software.
    # OSI Certified is a certification mark of the Open Source Initiative.
    #
    # Copyright (c) 2006, Enthought, Inc.
    # All rights reserved.
    #
    # Redistribution and use in source and binary forms, with or without
    # modification, are permitted provided that the following conditions are met:
    #
    # * Redistributions of source code must retain the above copyright notice, this
    #   list of conditions and the following disclaimer.
    # * Redistributions in binary form must reproduce the above copyright notice,
    #   this list of conditions and the following disclaimer in the documentation
    #   and/or other materials provided with the distribution.
    # * Neither the name of Enthought, Inc. nor the names of its contributors may
    #   be used to endorse or promote products derived from this software without
    #   specific prior written permission.
    #
    # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
    # ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    # WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    # DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
    # ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    # (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    # LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
    # ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    # (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    # SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    def _quat2transform(self):
        """
        Transform a unit quaternion into its corresponding rotation matrix (to
        be applied on the right side).

        :returns: Nx3x3 transform matrix
        :rtype: numpy array

        """
        q = self.q
        if q.ndim == 1:
            q = q[np.newaxis]

        x, y, z, w = q[:, 0], q[:, 1], q[:, 2], q[:, 3]
        xx2 = 2 * x * x
        yy2 = 2 * y * y
        zz2 = 2 * z * z
        xy2 = 2 * x * y
        wz2 = 2 * w * z
        zx2 = 2 * z * x
        wy2 = 2 * w * y
        yz2 = 2 * y * z
        wx2 = 2 * w * x

        rmat = np.empty((len(q), 3, 3), float)
        rmat[:, 0, 0] = 1.0 - yy2 - zz2
        rmat[:, 0, 1] = xy2 - wz2
        rmat[:, 0, 2] = zx2 + wy2
        rmat[:, 1, 0] = xy2 + wz2
        rmat[:, 1, 1] = 1.0 - xx2 - zz2
        rmat[:, 1, 2] = yz2 - wx2
        rmat[:, 2, 0] = zx2 - wy2
        rmat[:, 2, 1] = yz2 + wx2
        rmat[:, 2, 2] = 1.0 - xx2 - yy2

        return rmat

    def _equatorial2quat(self):
        """Dummy method to return return quat.

        :returns: quaternion
        :rtype: Quat

        """
        return self._transform2quat()

    def _equatorial2transform(self):
        """Construct the transform/rotation matrix from RA,Dec,Roll

        :returns: transform matrix
        :rtype: Nx3x3 numpy array

        """
        ra = np.radians(self._get_ra())
        dec = np.radians(self._get_dec())
        roll = np.radians(self._get_roll())
        ca = np.cos(ra)
        sa = np.sin(ra)
        cd = np.cos(dec)
        sd = np.sin(dec)
        cr = np.cos(roll)
        sr = np.sin(roll)
        # This is the transpose of the transformation matrix (related to
        # translation of original perl code
        rmat = np.array(
            [
                [ca * cd, sa * cd, sd],
                [-ca * sd * sr - sa * cr, -sa * sd * sr + ca * cr, cd * sr],
                [-ca * sd * cr + sa * sr, -sa * sd * cr - ca * sr, cd * cr],
            ]
        )

        return rmat.transpose()

    def _transform2quat(self):
        """Construct quaternions from the transform/rotation matrices

        :returns: quaternions formed from transform matrices
        :rtype: numpy array
        """

        transform = self.transform
        if transform.ndim == 2:
            transform = transform[np.newaxis]
        T = transform.transpose(0, 2, 1)
        # Code was copied from perl PDL code that uses backwards index ordering
        den = np.array(
            [
                1.0 + T[:, 0, 0] - T[:, 1, 1] - T[:, 2, 2],
                1.0 - T[:, 0, 0] + T[:, 1, 1] - T[:, 2, 2],
                1.0 - T[:, 0, 0] - T[:, 1, 1] + T[:, 2, 2],
                1.0 + T[:, 0, 0] + T[:, 1, 1] + T[:, 2, 2],
            ]
        )

        half_rt_q_max = 0.5 * np.sqrt(np.max(den, axis=0))
        max_idx = np.argmax(den, axis=0)
        poss_quat = np.zeros((4, len(T), 4))
        denom = 4.0 * half_rt_q_max
        poss_quat[0] = np.transpose(
            np.array(
                [
                    half_rt_q_max,
                    (T[:, 1, 0] + T[:, 0, 1]) / denom,
                    (T[:, 2, 0] + T[:, 0, 2]) / denom,
                    -(T[:, 2, 1] - T[:, 1, 2]) / denom,
                ]
            )
        )
        poss_quat[1] = np.transpose(
            np.array(
                [
                    (T[:, 1, 0] + T[:, 0, 1]) / denom,
                    half_rt_q_max,
                    (T[:, 2, 1] + T[:, 1, 2]) / denom,
                    -(T[:, 0, 2] - T[:, 2, 0]) / denom,
                ]
            )
        )
        poss_quat[2] = np.transpose(
            np.array(
                [
                    (T[:, 2, 0] + T[:, 0, 2]) / denom,
                    (T[:, 2, 1] + T[:, 1, 2]) / denom,
                    half_rt_q_max,
                    -(T[:, 1, 0] - T[:, 0, 1]) / denom,
                ]
            )
        )
        poss_quat[3] = np.transpose(
            np.array(
                [
                    -(T[:, 2, 1] - T[:, 1, 2]) / denom,
                    -(T[:, 0, 2] - T[:, 2, 0]) / denom,
                    -(T[:, 1, 0] - T[:, 0, 1]) / denom,
                    half_rt_q_max,
                ]
            )
        )

        q = np.zeros((len(T), 4))
        for idx in range(0, 4):
            max_match = max_idx == idx
            q[max_match] = poss_quat[idx][max_match]

        return q

    def __div__(self, quat2):
        """
        Divide one quaternion by another (or divide N quaternions by N quaternions)

        Example usage::

         >>> q1 = Quat((20,30,40))
         >>> q2 = Quat((30,40,50))
         >>> q = q1 / q2

        Performs the operation as q1 * inverse(q2) which is equivalent to
        the inverse(q2) transform followed by the q1 transform.  See the
        __mul__ operator help for more explanation on composing quaternions.

        :returns: product q1 * inverse q2
        :rtype: Quat

        """
        return self * quat2.inv()

    __truediv__ = __div__

    def __mul__(self, quat2):
        """
        Multiply quaternion by another.

        Quaternion composition as a multiplication q = q1 * q2 is equivalent to
        applying the q2 transform followed by the q1 transform.  Another way to
        express this is::

          q = Quat(numpy.dot(q1.transform, q2.transform))

        (though numpy.dot is not used because it is awkward in the vector case
        when the transforms are of the shape Nx3x3)

        Example usage::

          >>> q1 = Quat((20,30,0))
          >>> q2 = Quat((0,0,40))
          >>> (q1 * q2).equatorial
          array([20., 30., 40.])

        This example first rolls about X by 40 degrees, then rotates that rolled frame
        to RA=20 and Dec=30.  Doing the composition in the other order does a roll about
        (the original) X-axis of the (RA, Dec) = (20, 30) frame, yielding a non-intuitive
        though correct result::

          >>> (q2 * q1).equatorial
          array([ 353.37684725,   34.98868888,   47.499696  ])

        :returns: product q1 * q2
        :rtype: Quat

        """
        q1 = self.q
        if q1.ndim == 1:
            q1 = q1[np.newaxis]
        q2 = quat2.q
        if q2.ndim == 1:
            q2 = q2[np.newaxis]
        mult = np.zeros((len(q1), 4))
        mult[:, 0] = (
            q1[:, 3] * q2[:, 0]
            - q1[:, 2] * q2[:, 1]
            + q1[:, 1] * q2[:, 2]
            + q1[:, 0] * q2[:, 3]
        )
        mult[:, 1] = (
            q1[:, 2] * q2[:, 0]
            + q1[:, 3] * q2[:, 1]
            - q1[:, 0] * q2[:, 2]
            + q1[:, 1] * q2[:, 3]
        )
        mult[:, 2] = (
            -q1[:, 1] * q2[:, 0]
            + q1[:, 0] * q2[:, 1]
            + q1[:, 3] * q2[:, 2]
            + q1[:, 2] * q2[:, 3]
        )
        mult[:, 3] = (
            -q1[:, 0] * q2[:, 0]
            - q1[:, 1] * q2[:, 1]
            - q1[:, 2] * q2[:, 2]
            + q1[:, 3] * q2[:, 3]
        )
        return Quat(mult)

    def inv(self):
        """
        Invert the quaternion

        :returns: inverted quaternion
        :rtype: Quat
        """
        q = self.q
        if q.ndim == 1:
            q = q[np.newaxis]
        return Quat(np.array([q[:, 0], q[:, 1], q[:, 2], -1.0 * q[:, 3]]).transpose())


def normalize(array):
    """
    Normalize a 4 (or Nx4) element array/list/numpy.array for use as a quaternion

    :param quat_array: 4 or Nx4 element list/array
    :returns: normalized array
    :rtype: numpy array

    """
    quat = np.array(array)
    if quat.ndim == 1:
        return quat / np.sqrt(np.dot(quat, quat))
    elif quat.ndim == 2:
        return quat / np.sqrt(np.sum(quat * quat, axis=-1)[:, np.newaxis])
    else:
        raise TypeError("Input must be 1 or 2d")
