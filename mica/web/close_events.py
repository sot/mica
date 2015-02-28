from kadi import events
from Chandra.Time import DateTime

# 30 minute pad for these events
interest_pad = 30 * 60;
interest = [events.normal_suns,
            events.safe_suns,
            events.scs107s,
            events.caps]

def close_events(start, stop):
    es = []
    for etype in interest:
        e = etype.filter(start=DateTime(start).secs - interest_pad,
                         stop=DateTime(stop).secs + interest_pad)
        if e.count():
            es.append(e)
    me = events.major_events.filter(DateTime(start) - 1,
                                    DateTime(stop) + 1)
    if me.count():
        es.extend(me)
    return es


