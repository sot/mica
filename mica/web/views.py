from django.views.generic import TemplateView

from kadi import events
from kadi.events.views import BaseView
from Chandra.Time import DateTime

class IndexView(BaseView, TemplateView):
    template_name = 'mica/index.html'

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super(IndexView, self).get_context_data(**kwargs)

        obsid = self.request.GET.get('obsid_or_date', None)
        if obsid is not None:
            try:
                obsid = int(obsid)
            except:
                try:
                    obsids = events.obsids.filter(start=obsid)
                    obsid = obsids[0].obsid
                except:
                    obsid = None

        context['obsid'] = obsid or ''

        if obsid:
            obsid = format(obsid, '05d')
            url = ('https://icxc.harvard.edu/aspect/mica_reports/{}/{}/index.html'
                   .format(obsid[:2], obsid))
            context['mica_url'] = url

        return context


class AcqView(BaseView, TemplateView):
    template_name = 'mica/acq.html'

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super(AcqView, self).get_context_data(**kwargs)

        obsid = self.request.GET.get('obsid', None)
        if obsid is not None:
            try:
                obsid = int(obsid)
            except:
                obsid = None

        context['obsid'] = obsid or ''

        if obsid:
            import mica.web.pcad_table
            obsid = format(obsid, '05d')
            pcad_data = mica.web.pcad_table.get_acq_table(obsid)
            context['pcad_data'] = pcad_data

        return context


class TelemTableView(BaseView, TemplateView):
    template_name = 'mica/telem_table.html'

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super(TelemTableView, self).get_context_data(**kwargs)
        context['errors'] = []

        start_time = self.request.GET.get('start_time', None)
        stop_time = self.request.GET.get('stop_time', None)
        if start_time is not None:
            try:
                start_time = DateTime(start_time)
                start_time.date
            except:
                start_time = None
                context['start_time_date'] = ''
        if stop_time is not None:
            try:
                stop_time = DateTime(stop_time)
                stop_time.date
            except:
                stop_time = None
                context['stop_time_date'] = ''

        if start_time and stop_time is None:
            stop_time = DateTime(DateTime(start_time).secs + 25)
            start_time = DateTime(DateTime(start_time).secs - 25)


        if start_time and stop_time:
            context['start_time_date'] = start_time.date
            context['stop_time_date'] = stop_time.date
            import mica.web.close_events
            close_events = mica.web.close_events.close_events(start_time,
                                                              stop_time)
            print close_events
            context['events'] = close_events
            import mica.web.telem_table
            pcad_data = mica.web.telem_table.get_telem_table(start_time, stop_time)
            context['pcad_data'] = pcad_data

        return context
