# Licensed under a 3-clause BSD style license - see LICENSE.rst
from django.views.generic import TemplateView
from kadi import events
from kadi.events.views import BaseView


class IndexView(BaseView, TemplateView):
    template_name = "mica/index.html"

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super(IndexView, self).get_context_data(**kwargs)

        obsid = self.request.GET.get("obsid_or_date", None)
        if obsid is not None:
            try:
                obsid = int(obsid)
            except:
                try:
                    obsids = events.obsids.filter(start=obsid)
                    obsid = obsids[0].obsid
                except:
                    obsid = None

        context["obsid"] = obsid or ""

        if obsid:
            obsid = format(obsid, "05d")
            url = (
                "https://icxc.harvard.edu/aspect/mica_reports/{}/{}/index.html".format(
                    obsid[:2], obsid
                )
            )
            context["mica_url"] = url

        return context


class StarHistView(BaseView, TemplateView):
    template_name = "mica/star_hist.html"

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super(StarHistView, self).get_context_data(**kwargs)

        agasc_id = self.request.GET.get("agasc_id", None)
        if agasc_id is not None:
            try:
                agasc_id = int(agasc_id)
            except:
                agasc_id = None
        start = self.request.GET.get("start", None)
        stop = self.request.GET.get("stop", None)
        if start == "":
            start = None
        if stop == "":
            stop = None

        context["agasc_id"] = agasc_id or ""
        context["start"] = start or ""
        context["stop"] = stop or ""

        if agasc_id:
            import agasc
            from agasc.agasc import IdNotFound

            import mica.web.star_hist

            try:
                agasc_info = agasc.get_star(agasc_id, agasc_file="miniagasc_*")
                context["star_info"] = [
                    (key, agasc_info[key]) for key in agasc_info.dtype.names
                ]
            except IdNotFound:
                context["star_info"] = []
                pass
            acq_table, gui_table = mica.web.star_hist.get_star_stats(
                agasc_id, start, stop
            )
            if len(acq_table):
                context["acq_table"] = acq_table
            if len(gui_table):
                context["gui_table"] = gui_table
                reports_url = (
                    "https://cxc.cfa.harvard.edu/mta/ASPECT/agasc/supplement_reports/stars/"
                    + f"{int(agasc_id//1e7):03d}/{agasc_id}/index.html"
                )
                context["reports_url"] = reports_url
        return context


class AcqView(BaseView, TemplateView):
    template_name = "mica/acq.html"

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super(AcqView, self).get_context_data(**kwargs)

        obsid = self.request.GET.get("obsid", None)
        if obsid is not None:
            try:
                obsid = int(obsid)
            except:
                obsid = None

        context["obsid"] = obsid or ""

        if obsid:
            import mica.web.pcad_table

            obsid = format(obsid, "05d")
            pcad_data = mica.web.pcad_table.get_acq_table(obsid)
            context["pcad_data"] = pcad_data

        return context
