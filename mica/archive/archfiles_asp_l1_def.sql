CREATE TABLE archfiles (
  filename        text not null,
  filetime        int  not null,
  year            int not null,
  doy             int not null,
  tstart          float not null,
  tstop           float not null,
  checksum        text ,
  ascdsver        text ,
  caldbver        text ,
  content         text ,
  revision        int ,
  obsid           int not null,
  date            text not null,

  CONSTRAINT pk_archfiles PRIMARY KEY (filename)
);

CREATE INDEX idx_archfiles_filetime ON archfiles (filetime);
