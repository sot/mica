CREATE TABLE aspect_1_proc (
  aspect_1_id int not null,
  obsid  int not null,
  obi    int not null,
  revision int,
  ap_date text,
  isdefault int,
  obspar_version int,
  vv_complete int not null default 0,
  CONSTRAINT pk_archfiles PRIMARY KEY (aspect_1_id)
);

