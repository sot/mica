import os
import jinja2
from mica.stats.mag_stats import mag_stats
from astropy import table
import numpy as np
from cxotime import CxoTime
import matplotlib.pyplot as plt
import pandas as pd

OBS_STATS = None
AGASC_STATS = None


def _init(agasc_stats, obs_stats):
    global OBS_STATS
    global AGASC_STATS
    OBS_STATS = obs_stats
    AGASC_STATS = agasc_stats


def plot_agasc_id_single(agasc_stats, obs_stats, agasc_id, obsid=None, telem=None,
                  highlight_obsid=[], highlight_outliers=True, only_ok=True, draw_agasc_mag_stats=False, draw_obsid_mag_stats=False, draw_agasc_mag=False,
                  fit_ylim=True, title=None, draw_legend=False, ax=None):
    _init(agasc_stats, obs_stats)
    if title is not None:
        ax.set_title(title)
    if type(highlight_obsid) is not list and np.isscalar(highlight_obsid):
        highlight_obsid = [highlight_obsid]

    agasc_stat = AGASC_STATS[AGASC_STATS['agasc_id'] == agasc_id][0]
    obs_stats = OBS_STATS[OBS_STATS['agasc_id'] == agasc_id]
    if obsid:
        obs_stats = obs_stats[obs_stats['obsid'] == obsid]

    previous_axes = plt.gca()
    if ax is not None:
        plt.sca(ax)
    if ax is None:
        ax = plt.gca()
    if telem is None:
        telem = mag_stats.get_telemetry_by_agasc_id(agasc_id, ignore_exceptions=True)
        telem = mag_stats.add_obsid_info(telem, obs_stats)

    obsids = [obsid] if obsid else np.unique(telem['obsid'])

    timeline = telem[['times', 'mags', 'obsid', 'obsid_outlier']].copy()
    timeline['index'] = np.arange(len(timeline))
    timeline['mean'] = np.nan
    timeline['std'] = np.nan
    timeline['mag_mean'] = np.nan
    timeline['mag_std'] = np.nan
    for i, obsid in enumerate(np.unique(timeline['obsid'])):
        sel = (obs_stats['obsid'] == obsid)
        if draw_obsid_mag_stats and np.any(sel):
            timeline['mag_mean'][timeline['obsid'] == obsid] = obs_stats[sel]['mean'][0]
            timeline['mag_std'][timeline['obsid'] == obsid] = obs_stats[sel]['std'][0]
    timeline = timeline.as_array()

    ok = (telem['AOACASEQ'] == 'KALM') & (telem['AOACIIR'] == 'OK') & (telem['AOACISP'] == 'OK') & (telem['AOPCADMD'] == 'NPNT')
    ok = ok & (telem['dr'] < 3)

    highlighted = np.zeros(len(timeline['times']), dtype=bool)
    if highlight_obsid:
        highlighted = highlighted | np.in1d(timeline['obsid'], highlight_obsid)
    if highlight_outliers:
        highlighted = highlighted | timeline['obsid_outlier']

    limits = {}
    for i, obsid in enumerate(obsids):
        limits[obsid] = (timeline['index'][timeline['obsid'] == obsid].min(),
                         timeline['index'][timeline['obsid'] == obsid].max())

        if not only_ok and np.any((timeline['obsid'] == obsid) & ~ok):
            plt.scatter(timeline['index'][(timeline['obsid'] == obsid) & ~ok],
                        timeline[(timeline['obsid'] == obsid) & ~ok]['mags'],
                        s=10, marker='.', color='r')
        if np.any((timeline['obsid'] == obsid) & ok & highlighted):
            plt.scatter(timeline['index'][(timeline['obsid'] == obsid) & ok & highlighted],
                        timeline[(timeline['obsid'] == obsid) & ok & highlighted]['mags'],
                        s=10, marker='.', color='orange')
        if np.any((timeline['obsid'] == obsid) & ok & ~highlighted):
            plt.scatter(timeline['index'][(timeline['obsid'] == obsid) & ok & ~highlighted],
                        timeline[(timeline['obsid'] == obsid) & ok & ~highlighted]['mags'],
                        s=10, marker='.', color='k')
        sel = (obs_stats['obsid'] == obsid)
        if draw_obsid_mag_stats and np.sum(sel):
            label = '' if i else 'mag$_{OBSID}$'
            o = (timeline['obsid'] == obsid)
            if np.isfinite(timeline['mag_mean'][o][0]) and np.isfinite(timeline['mag_std'][o][0]):
                mag_mean_minus = timeline['mag_mean'][o][0] - timeline['mag_std'][o][0]
                mag_mean_plus = timeline['mag_mean'][o][0] + timeline['mag_std'][o][0]
                ax.plot(timeline[o]['index'],
                        timeline[o]['mag_mean'],
                        linewidth=2, color='orange', label=label)
                ax.fill_between([timeline['index'][o][0], timeline['index'][o][-1]],
                                [mag_mean_minus, mag_mean_minus],
                                [mag_mean_plus, mag_mean_plus],
                                color='orange', alpha=0.1, zorder=100)
            else:
                ax.plot([], [], linewidth=2, color='orange', label=label)
    if fit_ylim and agasc_stat['t_std_dr3'] > 0:
        ax.set_ylim((agasc_stat['t_mean_dr3'] - 6 * agasc_stat['t_std_dr3'],
                     agasc_stat['t_mean_dr3'] + 6 * agasc_stat['t_std_dr3']))
    
    sorted_obsids = sorted(limits.keys(), key=lambda l: limits[l][1])
    for i, obsid in enumerate(sorted_obsids):
        (tmin, tmax) = limits[obsid]
        ax.plot([tmin, tmin], ax.get_ylim(), ':', color='purple', scaley=False)
        shift = 0.07*(ax.get_ylim()[1] - ax.get_ylim()[0])*(1 + i%3)
        ax.text(np.mean([tmin, tmax]), ax.get_ylim()[0] + shift, f'{obsid}',
                verticalalignment='top', horizontalalignment='center')
    if limits:
        tmax = max([v[1] for v in limits.values()])
        ax.plot([tmax, tmax], ax.get_ylim(), ':', color='purple', scaley=False)
    
    ylim = ax.get_ylim()
    xlim = ax.get_xlim()
    ax.set_xlim(xlim)

    if draw_agasc_mag:
        mag_aca = np.mean(agasc_stat['mag_aca'])
        ax.plot(xlim, [mag_aca, mag_aca], label='mag$_{AGASC}$',
                color='green', scalex=False, scaley=False)

    if (draw_agasc_mag_stats and
            np.isfinite(agasc_stat['t_mean_dr3']) and agasc_stat['t_mean_dr3'] > 0):
        mag_weighted_mean = agasc_stat['t_mean_dr3']
        mag_weighted_std = agasc_stat['t_std_dr3']
        ax.plot(ax.get_xlim(), [mag_weighted_mean, mag_weighted_mean],
                label='mag', color='r', scalex=False)
        ax.fill_between(xlim,
                        [mag_weighted_mean - mag_weighted_std, mag_weighted_mean - mag_weighted_std],
                        [mag_weighted_mean + mag_weighted_std, mag_weighted_mean + mag_weighted_std],
                        color='r', alpha=0.1)
    
    if draw_legend:
        ax.set_xlim((xlim[0], xlim[1] + 0.1*(xlim[1] - xlim[0])))
        plt.legend(loc='center right')
    plt.sca(previous_axes)


def plot_flags(telemetry, ax=None, obsid=None):

    if ax is None:
        ax = plt.gca()

    timeline = telemetry[['times', 'mags', 'obsid', 'obsid_ok', 'dr', 'AOACFCT',
                          'AOACASEQ', 'AOACIIR', 'AOACISP', 'AOPCADMD',
                          ]]
    timeline['x'] = np.arange(len(timeline))
    timeline['y'] = np.ones(len(timeline))
    timeline = timeline.as_array()

    if obsid:
        timeline = timeline[timeline['obsid'] == obsid]

    obsids = np.unique(timeline['obsid'])

    limits = {}
    for i, obsid in enumerate(obsids):
        limits[obsid] = (timeline['x'][timeline['obsid'] == obsid].min(),
                         timeline['x'][timeline['obsid'] == obsid].max())

    all_ok = ((timeline['AOACASEQ'] == 'KALM') &
              (timeline['AOPCADMD'] == 'NPNT') &
              (timeline['AOACFCT'] == 'TRAK') &
              (timeline['AOACIIR'] == 'OK') &
              (timeline['AOACISP'] == 'OK') &
              (timeline['dr'] < 3) &
              timeline['obsid_ok']
              )
    flags = [
        ('OBS not OK', ~timeline['obsid_ok']),
        ('dr > 3', ((timeline['AOACASEQ'] == 'KALM') &
                    (timeline['AOPCADMD'] == 'NPNT') &
                    (timeline['AOACFCT'] == 'TRAK') &
                    (timeline['dr'] >= 3))),
        ('Ion. rad.', (timeline['AOACIIR'] != 'OK')),
        ('Sat. pixel.', (timeline['AOACISP'] != 'OK')),
        ('not track', ((timeline['AOACASEQ'] == 'KALM') &
                       (timeline['AOPCADMD'] == 'NPNT') &
                       (timeline['AOACFCT'] != 'TRAK'))),
        ('not Kalman', ((timeline['AOACASEQ'] != 'KALM') |
                        (timeline['AOPCADMD'] != 'NPNT'))),
    ]


    ok = [f[1] for f in flags]
    labels = [f[0] for f in flags]
    ticks = [i for i in range(len(flags))]

    for i in range(len(ok)):
        ax.scatter(timeline['x'][ok[i]], ticks[i]*timeline['y'][ok[i]], s=4, marker='.', color='k')
    ax.set_yticklabels(labels)
    ax.set_yticks(ticks)
    ax.set_ylim((-1, ticks[-1]+1))
    ax.grid(True, axis='y', linestyle=':')

    sorted_obsids = sorted(limits.keys(), key=lambda l: limits[l][1])
    for i, obsid in enumerate(sorted_obsids):
        (tmin, tmax) = limits[obsid]
        ax.axvline(tmin, linestyle=':', color='purple')
    if limits:
        ax.axvline(tmax, linestyle=':', color='purple')

    dr = timeline['dr'].copy()
    dr[dr > 10] = 10
    from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
    divider = make_axes_locatable(ax)
    ax_dr = divider.append_axes("bottom", size='25%', pad=0., sharex=ax)
    ax_dr.set_ylabel('dr')
    ax_dr.scatter(timeline['x'][all_ok & (dr < 10)], dr[all_ok & (dr < 10)], s=3, marker='.', color='k')
    ax_dr.scatter(timeline['x'][all_ok & (dr >= 10)], dr[all_ok & (dr >= 10)], s=3, marker='^', color='k')
    ax_dr.scatter(timeline['x'][~all_ok & (dr < 10)], dr[~all_ok & (dr < 10)], s=3, marker='.', color='r')
    ax_dr.scatter(timeline['x'][~all_ok & (dr >= 10)], dr[~all_ok & (dr >= 10)], s=3, marker='^', color='r')
    ax_dr.set_ylim((-0.5, 10.5))
    ax_dr.set_yticks([0., 2.5, 5, 7.5, 10], minor=True)
    ax_dr.set_yticks([0., 5, 10], minor=False)
    ax_dr.grid(True, axis='y', linestyle=':')

    for i, obsid in enumerate(sorted_obsids):
        (tmin, tmax) = limits[obsid]
        ax_dr.axvline(tmin, linestyle=':', color='purple')
    if limits:
        ax_dr.axvline(tmax, linestyle=':', color='purple')


def plot_set(agasc_stats, obs_stats, agasc_id, args, telem=None, filename=None):
    if not args:
        return
    if telem is None:
        telem = mag_stats.get_telemetry_by_agasc_id(agasc_id, ignore_exceptions=True)
        telem = mag_stats.add_obsid_info(telem, obs_stats[obs_stats['agasc_id'] == agasc_id])
    fig, ax = plt.subplots(len(args), 1, figsize=(15, 3.5*len(args)))
    if len(args) == 1:
        ax = [ax]
    ax[0].set_title(f'AGASC ID {agasc_id}')

    for i, kwargs in enumerate(args):
        if 'type' in kwargs and kwargs['type'] == 'flags':
            o = kwargs['obsid'] if 'obsid' in kwargs else None
            plot_flags(telem, ax[i], obsid=o)
            if i:
                ax[i].set_xlim(ax[i-1].get_xlim())
        else:
            plot_agasc_id_single(agasc_stats, obs_stats, agasc_id, telem=telem, ax=ax[i], **kwargs)

    plt.tight_layout()
    if filename is not None:
        fig.savefig(filename)

    return fig


STAR_REPORT_BOOTSTRAP = """<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
    <link rel="stylesheet"
          href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css"
          integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm"
          crossorigin="anonymous">
  </head>
  <body>

    <div class="container">
      <h1> AGASC ID {{ agasc_stats.agasc_id }} </h1>
      <h3> Info </h3>
      <div class="row">
        <div class="col-md">
          <table class="table table-bordered table-sm">
            <tr> <td style="width: 30%"> Last Obs. </td> <td style="width: 30%"> {{ agasc_stats.last_obs }} </td> </tr>
            <tr>
              <td style="width: 30%"> mag<sub>catalog</sub> </td>
              <td style="width: 30%"> {{ "%.2f" | format(agasc_stats.mag_aca) }} &#177; {{ "%.2f" | format(agasc_stats.mag_aca_err) }} </td>
            </tr>
            <tr>
              <td> mag<sub>3 arcsec </sub> </td>
              <td> {{ "%.2f" | format(agasc_stats.t_mean_dr3) }} &#177; {{ "%.2f" | format(agasc_stats.t_std_dr3) }} </td>
            </tr>
            <tr>
              <td> mag<sub>5 arcsec </sub> </td>
              <td> {{ "%.2f" | format(agasc_stats.t_mean_dr5) }} &#177; {{ "%.2f" | format(agasc_stats.t_std_dr5) }} </td>
            </tr>
          </table>
        </div>
        <div class="col-md">
          <table class="table table-bordered table-sm">
            <tr>
              <td> N<sub>obs</sub> </td>
              <td> {{ agasc_stats.n_obsids }} <span{%- if agasc_stats.n_obs_bad %} style="color:red;"{% endif -%}> ({{ agasc_stats.n_obs_bad }} bad) <span> </td>
            </tr>
            <tr> <td> f<sub>ok</sub> </td> <td> {{ "%.1f" | format(100*agasc_stats.f_ok) }}%  </td> </tr>
            <tr> <td> f<sub>3 arcsec</sub> </td> <td> {{ "%.1f" | format(100*agasc_stats.f_dr3) }}% </td> </tr>
            <tr> <td> f<sub>5 arcsec</sub> </td> <td> {{ "%.1f" | format(100*agasc_stats.f_dr5) }}% </td> </tr>
          </table>
        </div>
      </div>


      <h3> Timeline </h3>
      <img src="mag_stats.png" width="100%"/>

      <h3> Observation Info </h3>
      <table  class="table table-hover">
        <tr>
          <th data-toggle="tooltip" data-placement="top" title="OBSID"> OBSID </th>
          <th data-toggle="tooltip" data-placement="top" title="MP starcat time"> Time </th>
          <th data-toggle="tooltip" data-placement="top" title="Pixel row"> Row </th>
          <th data-toggle="tooltip" data-placement="top" title="Pixel column"> Col </th>
          <!-- th data-toggle="tooltip" data-placement="top" data-html="true" title="Observation is considered in the calculation <br/> n &gt; 10 <br/>f_ok &gt; 0.3 <br/> &langle; &delta; <sub>mag</sub> &rangle; <sub>100s</sub>  < 1"> OK </th -->
          <th data-toggle="tooltip" data-placement="top" data-html="true" title="Number of time samples"> N </th>
          <th data-toggle="tooltip" data-placement="top" data-html="true" title="Number of time samples considered as 'tracking' <br/> AOACASEQ == 'KALM' <br/> AOACIIR == 'OK' <br/> AOACISP == 'OK' <br/> AOPCADMD == 'NPNT' <br/> AOACFCT == 'TRAK' <br/> OBS_OK)"> N<sub>ok</sub> </th>
          <th data-toggle="tooltip" data-placement="top" data-html="true" title="Number of outlying samples"> N<sub>out</sub> </th>
          <th data-toggle="tooltip" data-placement="top" data-html="true" title="Tracking time as a fraction of the total time <br/> AOACASEQ == 'KALM' <br/> AOACIIR == 'OK' <br/> AOACISP == 'OK' <br/> AOPCADMD == 'NPNT' <br/> AOACFCT == 'TRAK' <br/> OBS_OK)"> f<sub>track</sub> </th>
          <th data-toggle="tooltip" data-placement="top" data-html="true" title="Fraction of tracking time within 3 arcsec of target"> f<sub>dr3</sub> </th>
          <th data-toggle="tooltip" data-placement="top" data-html="true" title="Fraction of tracking time within 5 arcsec of target"> f<sub>dr5</sub> </th>
          <th data-toggle="tooltip" data-placement="top" data-html="true" title="Time where slot is tracking and target within 3 arcsec <br/> as fraction of total time"> f<sub>ok</sub> </th>
          <th data-toggle="tooltip" data-placement="top" data-html="true" title="100-second Rolling mean of mag - &langle; mag &rangle;"> &langle; &delta; <sub>mag</sub> &rangle; <sub>100s</sub>  </th>
          <th data-toggle="tooltip" data-placement="top" data-html="true" title="Mean magnitude"> &langle; mag &rangle; </th>
          <th data-toggle="tooltip" data-placement="top" data-html="true" title="Magnitude uncertainty"> &sigma;<sub>mag</sub> </th>
        </tr>
        {%- for s in obs_stats %}
        <tr {%- if not s.obsid_ok %} class="table-danger" {% endif %}>
          <td> <a href="https://web-kadi.cfa.harvard.edu/mica/?obsid_or_date={{ s.obsid }}"> {{ s.obsid }} </td>
          <td> {{ s.mp_starcat_time }} </td>
          <td> {{ "%.1f" | format(s.row) }} </td>
          <td> {{ "%.1f" | format(s.col) }} </td>
          <!-- td> {{ s.obsid_ok }} </td -->
          <td> {{ s.n }} </td>
          <td> {{ s.n_ok }} </td>
          <td> {{ s.outliers }} </td>
          <td> {{ "%.1f" | format(100*s.f_track) }}% </td>
          <td> {{ "%.1f" | format(100*s.f_dr3) }}% </td>
          <td> {{ "%.1f" | format(100*s.f_dr5) }}% </td>
          <td> {{ "%.1f" | format(100*s.f_ok) }}% </td>
          <td> {{ "%.2f" | format(s.lf_variability_100s) }} </td>
          <td> {{ "%.2f" | format(s.t_mean) }} </td>
          <td> {{ "%.2f" | format(s.t_mean_err) }} </td>
        </tr>
        {%- endfor %}
      </table>

    </div>
  </body>
</html>
"""

STAR_REPORT_TABULATOR = """<html>
  <head>
    <link href="{{ static_dir }}/tabulator/dist/css/semantic-ui/tabulator_semantic-ui.min.css"
          rel="stylesheet">
    <script type="text/javascript"
            src="{{ static_dir }}/tabulator/dist/js/tabulator.min.js"></script>
    <style type="text/css">
      body { padding: 30px; font-size: 20px }
      #stats_table {
        font-size: large;
        border: 1px solid #333;
        border-radius: 10px;
      }
      #obs_table {
        font-size: large;
        border: 1px solid #333;
        border-radius: 10px;
      }
    </style>
  </head>
  <body>

    <h1> Magnitude Statistics </h1>
    <h2> AGASC ID 31989744 </h2>

    <div id="stats_table" style="width:auto"></div>

    <img src="mag_stats.png"/>

    <div id="obs_table" style="width:auto"></div>

    <script type="text/javascript">
    var stats_table = new Tabulator("#stats_table", {
        responsiveLayout:"hide",
        layout:"fitDataFill",
        columns:[
          {
            field: "field",
          },
          {
            field: "value",
          }
        ],

      });

    var stats_data = [
      {
        field: "AGASC ID",
        value: "{{ agasc_stats.agasc_id }}"
      },
      {
        field: "last_obs",
        value: "{{ agasc_stats.last_obs }}"
      },
      {
        field: "n_obsids",
        value: "{{ agasc_stats.n_obsids }}"
      },
      {
        field: "n_obsids_ok",
        value: "{{ agasc_stats.n_obsids_ok }}"
      },
      {
        field: "f_ok",
        value: "{{ '%.1f' | format(100*agasc_stats.f_ok) }}%"
      },
      {
        field: "f_dr3",
        value: "{{ '%.1f' | format(100*agasc_stats.f_dr3) }}%"
      },
      {
        field: "f_dr5",
        value: "{{ '%.1f' | format(100*agasc_stats.f_dr5) }}%"
      },
      {
        field: "mag_aca",
        value: "{{ '%.2f' | format(agasc_stats.mag_aca) }}"
      },
      {
        field: "mag_aca_err",
        value: "{{ '%.2f' | format(agasc_stats.mag_aca_err) }}"
      },
      {
        field: "t_mean_dr3",
        value: "{{ '%.2f' | format(agasc_stats.t_mean_dr3) }}"
      },
      {
        field: "t_std_dr3",
        value: "{{ '%.2f' | format(agasc_stats.t_std_dr3) }}"
      },
      {
        field: "t_mean_dr5",
        value: "{{ '%.2f' | format(agasc_stats.t_mean_dr5) }}"
      },
      {
        field: "t_std_dr5",
        value: "{{ '%.2f' | format(agasc_stats.t_std_dr5) }}"
      },
    ];

    stats_table.setData(stats_data);

    var obs_table = new Tabulator("#obs_table", {
      dataTree:true,
      height: "auto",
      responsiveLayout:"hide",
      //layout:"fitDataFill",
      layout:"fitColumns",
      columns:[
        {field: "obsid", title: "OBSID",},
        {field: "mp_starcat_time", title: "time"},
        {field: "obsid_ok", title: "OK",},
        {field: "n", title: "N",},
        {field: "n_ok", title: "N<sub>ok</sub>",},
        {field: "outliers", title: "N<sub>out</sub>",},
        {field: "f_track", title: "f<sub>track</sub>",},
        {field: "f_ok", title: "f<sub>ok</sub>",},
        {field: "f_dr3", title: "f<sub>dr3</sub>",},
        {field: "f_dr5", title: "f<sub>dr5</sub>",},
        {field: "f_14", title: "f<sub>14</sub>",},
        {field: "lf_variability_100s", title: "&langle; &delta; <sub>mag</sub> &rangle; <sub>100s</sub> ",},
        {field: "t_mean", title: "&langle; mag &rangle;",},
        {field: "t_mean_err", title: "&sigma;<sub>mag</sub>",},
        ],
      });

      var obs_data = [
      {%- for s in obs_stats %}
        {
          obsid: "{{ s.obsid }}",
          mp_starcat_time: "{{ s.mp_starcat_time }}",
          obsid_ok: "{{ s.obsid_ok }}",
          n: "{{ s.n }}",
          n_ok: "{{ s.n_ok }}",
          outliers: "{{ s.outliers }}",
          f_track: "{{ "%.1f" | format(100*s.f_track) }}%",
          f_ok: "{{ "%.1f" | format(100*s.f_ok) }}%",
          f_dr3: "{{ "%.1f" | format(100*s.f_dr3) }}%",
          f_dr5: "{{ "%.1f" | format(100*s.f_dr5) }}%",
          f_14: "{{ "%.1f" | format(100*s.f_14) }}%",
          lf_variability_100s: "{{ "%.2f" | format(s.lf_variability_100s) }}",
          t_mean: "{{ "%.2f" | format(s.t_mean) }}",
          t_mean_err:         "{{ "%.2f" | format(s.t_mean_err) }}",
        },
      {%- endfor %}
      ];

     obs_table.setData(obs_data);
    </script>
  </body>
</html>
"""


def single_star_html_report(agasc_stats, obs_stats, agasc_id, directory='./mag_stats_reports',
                            static_dir='https://cxc.cfa.harvard.edu/mta/ASPECT/www_resources'):
    _init(agasc_stats, obs_stats)
    if np.sum(AGASC_STATS['agasc_id'] == agasc_id) == 0:
        return

    star_template = jinja2.Template(STAR_REPORT_BOOTSTRAP)

    if directory is None:
        directory = os.path.join(directory, 'stars', f'{agasc_id // 1e7:03.0f}', f'{agasc_id:.0f}')

    if not os.path.exists(directory):
        os.makedirs(directory)

    o = OBS_STATS[OBS_STATS['agasc_id'] == agasc_id]
    if len(o) == 0:
        raise Exception(f'agasc_id {agasc_id} has not observations')
    o.sort(keys=['mp_starcat_time'])
    s = AGASC_STATS[AGASC_STATS['agasc_id'] == agasc_id][0]
    s = {k: s[k] for k in s.colnames}
    s['n_obs_bad'] = s['n_obsids'] - s['n_obsids_ok']
    s['last_obs'] = ':'.join(o[-1]['mp_starcat_time'].split(':')[:4])

    # OBSIDs can be repeated
    obsids = list(np.unique(o[~o['obsid_ok']]['obsid']))

    args = [{'only_ok': False, 'draw_agasc_mag': True, 'draw_legend': True},
            {'title': 'mag stats',
             'only_ok': True,
             'highlight_obsid': obsids,
             'draw_obsid_mag_stats': True,
             'draw_agasc_mag_stats': True,
             'draw_legend': True
             },
            {'type': 'flags'}]
    for obsid in obsids:
        args.append({'obsid': obsid,
                     'only_ok': False,
                     'draw_obsid_mag_stats': True,
                     'draw_agasc_mag_stats': True,
                     'draw_legend': True
                     })
        args.append({'type': 'flags'})
    fig = plot_set(agasc_stats, obs_stats, agasc_id, args=args, filename=os.path.join(directory, f'mag_stats.png'))
    plt.close(fig)

    with open(os.path.join(directory, 'index.html'), 'w') as out:
        out.write(star_template.render(agasc_stats=s,
                                       obs_stats=o.as_array(),
                                       static_dir=static_dir))
    return os.path.join(directory, 'index.html')


RUN_REPORT_SIMPLE = """<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
    <link rel="stylesheet"
          href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css"
          integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm"
          crossorigin="anonymous">
  </head>
  <body>

    <div class="container">
    <h1> ACA Magnitude Statistics </h1>
    <h2> {{ info.report_date }} Update Report </h2>
    <table class="table table-sm">
      <tr>
        <td style="width: 50%"> Time range </td>
        <td style="width: 50%"> {{ info.tstart }} &ndash; {{ info.tstop }} </td>
      </tr>
      {%- for section in sections %}
      <tr>
        <td> <a href="#{{ section.id }}"> {{ section.title }} </a> </td>
        <td> {{ section.stars | length }} </td>
      </tr>
      {%- endfor %}    
      <tr>
        <td> {% if failures -%} <a href="#failures"> Failures </a>
             {%- else -%} Failures {%- endif %} </td>
        <td> {{ failures | length }} </td>
      </tr>
    </table>

    {%- for section in sections %}
    <a name="{{ section.id }}"> </a>
    <h3> {{ section.title }} </h3>
    <table class="table table-hover">
      <tr>
      <tr>
        <th data-toggle="tooltip" data-placement="top" title="ID in AGASC"> AGASC ID </th>
        <th data-toggle="tooltip" data-placement="top" title="Last time the star was observed"> Last Obs </th>
        <th data-toggle="tooltip" data-placement="top" title="Number of times the star has been observed"> n<sub>obs</sub> </th>
        <th data-toggle="tooltip" data-html="true" data-placement="top" title="Observations not included in calculation <br/> n &gt; 10 <br/>f_ok &gt; 0.3 <br/> &langle; &delta; <sub>mag</sub> &rangle; <sub>100s</sub>  < 1"> n<sub>bad</sub> </th>
        <th data-toggle="tooltip" data-html="true" data-placement="top" title="New observations not included in calculation <br/> n &gt; 10 <br/>f_ok &gt; 0.3 <br/> &langle; &delta; <sub>mag</sub> &rangle; <sub>100s</sub>  < 1"> n<sub>bad new</sub> </th>
        <th data-toggle="tooltip" data-placement="top" data-html="true" title="tracking time as fraction of total time: <br/> AOACASEQ == 'KALM' <br/> AOACIIR == 'OK' <br/> AOACISP == 'OK' <br/> AOPCADMD == 'NPNT' <br/> AOACFCT == 'TRAK' <br/> OBS_OK)"> f<sub>track</sub> </th>
        <th data-toggle="tooltip" data-placement="top" title="Fraction of the tracking time within 3 arcsec of target"> f<sub>3 arcsec</sub> </th>
        <th data-toggle="tooltip" data-placement="top" title="Fraction of the tracking time within 5 arcsec of target"> f<sub>5 arcsec</sub> </th>
        <th data-toggle="tooltip" data-placement="top" title="Magnitude in AGASC"> mag<sub>catalog</sub> </th>
        <th data-toggle="tooltip" data-placement="top" title="Magnitude observed"> mag<sub>obs</sub> </th>
        <th data-toggle="tooltip" data-placement="top" title="Difference between observed and catalog magnitudes"> &delta;<sub>mag cat</sub> </th>
        <th data-toggle="tooltip" data-placement="top" title="Difference between observed and catalog magnitudes, divided by catalog magnitude error"> &delta;<sub>mag</sub>/&sigma;<sub>mag</sub> </th>
        <th data-toggle="tooltip" data-placement="top" title="Variation in observed magnitude from the last version of AGASC supplement"> &delta;<sub>mag</sub> </th>
        <th data-toggle="tooltip" data-placement="top" title="Variation in observed magnitude standard deviation from the last version of AGASC supplement"> &delta;<sub>&sigma;</sub> </th>
        <th data-toggle="tooltip" data-placement="top" title="Color in AGASC"> color </th>
      </tr>
      {%- for star in section.stars %}
      <tr {%- if star.flag != '' %} class="table-{{ star.flag }}" {% endif %}>
        <td> {%- if star.agasc_id in star_reports %} <a href="{{ star_reports[star.agasc_id] }}/index.html"> {{ star.agasc_id }} </a> {% else %} {{ star.agasc_id }} {% endif %} </td>
        <td> {{ star.last_obs[:8] }} </td>
        <td> {{ star.n_obsids }}  </td>
        <td> {%- if star.n_obs_bad > 0 %} {{ star.n_obs_bad }} {% endif %} </td>
        <td> {%- if star.n_obs_bad > 0 %} {{ star.n_obs_bad_new }} {% endif %} </td>
        <td> {{ "%.1f" | format(100*star.f_ok) }}%  </td>
        <td> {{ "%.1f" | format(100*star.f_dr3) }}% </td>
        <td> {{ "%.1f" | format(100*star.f_dr5) }}% </td>
        <td {%- if star.selected_mag_aca_err %} class="table-info" {% endif %}> {{ "%.2f" | format(star.mag_aca) }} &#177; {{ "%.2f" | format(star.mag_aca_err) }}  </td>
        <td> {{ "%.2f" | format(star.mag_obs) }} &#177; {{ "%.2f" | format(star.mag_obs_err) }}  </td>
        <td {%- if star.selected_atol %} class="table-info" {% endif %}> {{ "%.2f" | format(star.delta) }}  </td>
        <td {%- if star.selected_rtol %} class="table-info" {% endif %}> {{ "%.2f" | format(star.sigma) }}  </td>
        <td> {% if star.new %} &ndash; {% else %}{{ "%.2f" | format(star.update_mag_aca) }}{% endif %}  </td>
        <td> {% if star.new %} &ndash; {% else %}{{ "%.2f" | format(star.update_mag_aca_err) }}{% endif %}  </td>
        <td {%- if star.selected_color %} class="table-info" {% endif %}> {{ "%.2f" | format(star.color) }}  </td>
      </tr>
      {%- endfor %}
    <table>
    {%- endfor %}

    <a name="failures"> </a>
    {%- if failures %}
    <h3> Failures </h3>
    <table class="table table-hover">
      <tr>
        <th> AGASC ID </th>
        <th> OBSID </th>
        <th> Message </th>
      </tr>
      {%- for failure in failures %}
      <tr>
        <td> {%- if failure.agasc_id in star_reports %} <a href="{{ star_reports[failure.agasc_id] }}/index.html"> {{ failure.agasc_id }} </a> {% else %} {{ failure.agasc_id }} {% endif %} </td>
        <td> {{ failure.obsid }} </td>
        <td> {{ failure.msg }} </td>
      </tr>
      {%- endfor %}
    </table>
    {% endif %}
    </div>
    
  <script src="https://code.jquery.com/jquery-3.2.1.slim.min.js" integrity="sha384-KJ3o2DKtIkvYIK3UENzmM7KCkRr/rE9/Qpg6aAZGJwFDMVNA/GpGFF93hXpG5KkN" crossorigin="anonymous"></script>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.9/umd/popper.min.js" integrity="sha384-ApNbgh9B+Y1QKtv3Rn7W3mgPxhU9K/ScQsAP7hUibX39j7fakFPskvXusvfa0b4Q" crossorigin="anonymous"></script>
  <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js" integrity="sha384-JZR6Spejh4U02d8jOt6vLEHfe/JQGiRRSQQxSfFWpi1MquVdAyjUar5+76PVCmYl" crossorigin="anonymous"></script>

  <script type="text/javascript">
    $(document).ready(function() {
    $("body").tooltip({ selector: '[data-toggle=tooltip]' });
});
  </script>
</html>
"""


def multi_star_html_report(agasc_stats, obs_stats, sections={}, new_stars=[], updated_stars=[], fails=[],
                           tstart=None, tstop=None, report_date=None,
                           directory='./mag_stats_reports', filename=None, include_all_stars=False,
                           make_single_reports=True):
    _init(agasc_stats, obs_stats)

    run_template = jinja2.Template(RUN_REPORT_SIMPLE)

    updated_star_ids = updated_stars['agasc_id'] if len(updated_stars) else []

    info = {
        'tstart': tstart if tstart else CxoTime(OBS_STATS['mp_starcat_time']).min().date,
        'tstop': tstop if tstop else CxoTime(OBS_STATS['mp_starcat_time']).max().date,
        'report_date': report_date if report_date else CxoTime.now().date
    }

    if filename is None:
        filename = f'mag_stats_{info["tstart"]}-{info["tstop"]}.html'

    # this is the list of agasc_id for which we will generate individual reports (if possible)
    agasc_ids = np.concatenate([np.array(s['stars'], dtype=int) for s in sections])
    if include_all_stars:
        sections.append({
            'id': 'other_stars',
            'title': 'Other Stars',
            'stars': AGASC_STATS['agasc_id'][~np.in1d(AGASC_STATS['agasc_id'], agasc_ids)]
        })
        agasc_ids = AGASC_STATS['agasc_id']
    failed_agasc_ids = [f['agasc_id'] for f in fails
                        if f['agasc_id'] and int(f['agasc_id']) in OBS_STATS['agasc_id']]
    agasc_ids = np.unique(np.concatenate([agasc_ids, failed_agasc_ids]))

    # this turns all None into '' in a new list of failures
    fails = [{k: '' if v is None else v for k, v in f.items()} for i, f in enumerate(fails)]

    # check how many observations were added in this run, and how many of those are ok
    new_obs = OBS_STATS[(OBS_STATS['mp_starcat_time'] >= info["tstart"]) &
                        (OBS_STATS['mp_starcat_time'] <= info["tstop"])]. \
        group_by('agasc_id')[['agasc_id', 'obsid', 'obsid_ok']]. \
        groups.aggregate(np.count_nonzero)[['agasc_id', 'obsid', 'obsid_ok']]
    new_obs['n_obs_bad_new'] = new_obs['obsid'] - new_obs['obsid_ok']

    # add some extra fields
    agasc_stats = AGASC_STATS.copy()
    if len(agasc_stats) == 0:
        return [], []
    all_agasc_ids = np.unique(np.concatenate([
        new_obs['agasc_id'],
        [f['agasc_id'] for f in fails]
    ]))

    assert np.all(np.in1d(agasc_stats['agasc_id'], all_agasc_ids)), 'Not all AGASC IDs are in new obs.'
    agasc_stats['n_obs_bad'] = agasc_stats['n_obsids'] - agasc_stats['n_obsids_ok']
    agasc_stats['flag'] = '          '
    if len(agasc_stats):
        agasc_stats = table.join(agasc_stats, new_obs[['agasc_id', 'n_obs_bad_new']],
                                 keys=['agasc_id'])
    agasc_stats['flag'][:] = ''
    agasc_stats['flag'][agasc_stats['n_obs_bad'] > 0] = 'warning'
    agasc_stats['flag'][agasc_stats['n_obs_bad_new'] > 0] = 'danger'
    agasc_stats['delta'] = (agasc_stats['t_mean_dr3'] - agasc_stats['mag_aca'])
    agasc_stats['sigma'] = (agasc_stats['t_mean_dr3'] - agasc_stats['mag_aca'])/agasc_stats['mag_aca_err']
    agasc_stats['new'] = True
    agasc_stats['new'][np.in1d(agasc_stats['agasc_id'], updated_star_ids)] = False
    agasc_stats['update_mag_aca'] = np.nan
    agasc_stats['update_mag_aca_err'] = np.nan
    agasc_stats['last_obs'] = CxoTime(agasc_stats['last_obs_time']).date
    if len(updated_stars):
        agasc_stats['update_mag_aca'][np.in1d(agasc_stats['agasc_id'], updated_star_ids)] = \
            updated_stars['mag_aca']
        agasc_stats['update_mag_aca_err'][np.in1d(agasc_stats['agasc_id'], updated_star_ids)] = \
            updated_stars['mag_aca_err']
    # make all individual star reports
    star_reports = {}
    for agasc_id in np.atleast_1d(agasc_ids):
        try:
            dirname = os.path.join(directory, 'stars', f'{agasc_id//1e7:03.0f}', f'{agasc_id:.0f}')
            if make_single_reports:
                single_star_html_report(agasc_stats, obs_stats, agasc_id, directory=dirname)
            if os.path.exists(dirname):
                star_reports[agasc_id] = dirname
        except mag_stats.MagStatsException:
            pass

    # remove empty sections, and set the star tables for each of the remaining sections
    sections = sections.copy()
    sections = [section for section in sections if len(section['stars'])]
    for section in sections:
        section['stars'] = agasc_stats[np.in1d(agasc_stats['agasc_id'],
                                               section['stars'])].as_array()

    # this is a hack
    star_reports = {i: os.path.relpath(star_reports[i], directory) for i in star_reports}
    # make report
    if not os.path.exists(directory):
        os.makedirs(directory)
    with open(os.path.join(directory, filename), 'w') as out:
        out.write(run_template.render(info=info,
                                      sections=sections,
                                      failures=fails,
                                      star_reports=star_reports))
