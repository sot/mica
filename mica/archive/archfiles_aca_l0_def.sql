CREATE TABLE archfiles (
  filename        text not null,
  filetime        int  not null,
  year            int not null,
  doy             int not null,
  tstart          float not null,
  tstop           float not null,
  startmjf        int ,
  startmnf        int ,
  stopmjf         int ,
  stopmnf         int ,
  checksum        text ,
  tlmver          text ,
  ascdsver        text ,
  revision        int ,
  date            text not null,
  rows            int not null,
  imgsize         int not null,
  slot            int not null,
  CONSTRAINT pk_archfiles PRIMARY KEY (filename)
);

CREATE INDEX idx_archfiles_filetime ON archfiles (filetime);
