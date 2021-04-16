#include "../gerda_ar39_pdf.hpp"

std::vector<int> channels = {
    0, 1, 2, 3, 4, 5, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
    23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 37, 38, 39, 40
};

void fccd_sweep() {

    double dlf = 0.2;

    for (auto ch : channels) {
        auto c = new TCanvas(Form("ch%i", ch), "c");
        c->cd();
        std::vector<TH1D*> histos;
        for (double fccd = 0.65; fccd <= 2.4; fccd += 0.05) {

            auto h = new TH1D(Form("channel %i, fccd = %.2f mm, dlf = %.2f", ch, fccd, dlf),
                    Form("channel %i, fccd = 0.65:0.05:2.4 mm, dlf = %.2f", ch, dlf), 400, 0, 400);

            for (int i = 1; i <= h->GetNbinsX(); ++i) {
                h->SetBinContent(i, gerda::ar39_pdf(ch, h->GetBinCenter(i), fccd, dlf));
            }
            h->Draw("c same");
            histos.push_back(h);
        }
        c->SaveAs(Form("plots/fccd/ch%i.pdf", ch));
        delete c;
        for (auto h : histos) delete h;
        histos.clear();
    }
}

void dlf_sweep() {

    double fccd = 1.;

    for (auto ch : channels) {
        auto c = new TCanvas(Form("ch%i", ch), "c");
        c->cd();
        std::vector<TH1D*> histos;
        for (double dlf = 1; dlf >= 0; dlf -= 0.05) {

            auto h = new TH1D(Form("channel %i, fccd = %.2f mm, dlf = %.2f", ch, fccd, dlf),
                    Form("channel %i, fccd = %.2f mm, dlf = 0:0.05:1", ch, fccd), 400, 0, 400);

            for (int i = 1; i <= h->GetNbinsX(); ++i) {
                h->SetBinContent(i, gerda::ar39_pdf(ch, h->GetBinCenter(i), fccd, dlf));
            }
            h->Draw("c same");
            histos.push_back(h);
        }
        c->SaveAs(Form("plots/dlf/ch%i.pdf", ch));
        delete c;
        for (auto h : histos) delete h;
        histos.clear();
    }
}

void test() {

    gROOT->SetBatch();
    fccd_sweep();
    dlf_sweep();
}
