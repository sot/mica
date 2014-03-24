CREATE TABLE report_proc (
  obsid  int not null,
  checked_date text not null,
  ocat_status text,
  report_status text,
  report_version text,
  vv_version text,
  vv_revision text,
  aspect_1_id text,
  long_term text,
  short_term text,
  starcheck text,
  CONSTRAINT pk_archfiles PRIMARY KEY (obsid, report_version)
);

