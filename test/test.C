#include "../include/gerda_ar39_pdf.hpp"

std::vector<int> channels = {
    0, 1, 2, 3, 4, 5, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
    23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 37, 38, 39, 40
};

std::map<int,std::string> det_names = {
    {0, "GD91A"},
    {1, "GD35B"},
    {2, "GD02B"},
    {3, "GD00B"},
    {4, "GD61A"},
    {5, "GD89B"},
    {6, "GD02D"},
    {7, "GD91C"},
    {8, "ANG5"},
    {9, "RG1"},
    {10, "ANG3"},
    {11, "GD02A"},
    {12, "GD32B"},
    {13, "GD32A"},
    {14, "GD32C"},
    {15, "GD89C"},
    {16, "GD61C"},
    {17, "GD76B"},
    {18, "GD00C"},
    {19, "GD35C"},
    {20, "GD76C"},
    {21, "GD89D"},
    {22, "GD00D"},
    {23, "GD79C"},
    {24, "GD35A"},
    {25, "GD91B"},
    {26, "GD61B"},
    {27, "ANG2"},
    {28, "RG2"},
    {29, "ANG4"},
    {30, "GD00A"},
    {31, "GD02C"},
    {32, "GD79B"},
    {33, "GD91D"},
    {34, "GD32D"},
    {35, "GD89A"},
    {36, "IC48B"},
    {37, "IC50B"},
    {38, "IC48A"},
    {39, "IC50A"},
    {40, "IC74A"}
};

auto freec = 10000;

void fccd_sweep(double dlf=0.2) {

    for (auto ch : channels) {
        auto c = new TCanvas(Form("ch%i", ch), "c");
        c->cd();
        std::vector<TH1D*> histos;

        // in grid
        for (double fccd = 0.65; fccd <= 2.4; fccd += 0.05) {

            auto h = new TH1D(Form("%s, fccd = %.2f mm, dlf = %.2f", det_names[ch].c_str(), fccd, dlf),
                    Form("%s, fccd = 0.65:0.025:2.4 mm, dlf = %.2f", det_names[ch].c_str(), dlf), 200, 0, 200);

            for (int i = 1; i <= h->GetNbinsX(); ++i) {
                h->SetBinContent(i, gerda::ar39_pdf(ch, h->GetBinCenter(i), fccd, dlf));
            }
            h->SetLineColor(kRed);
            h->Draw("l same");
            histos.push_back(h);
        }
        // out grid
        for (double fccd = 0.675; fccd <= 2.4; fccd += 0.05) {

            auto h = new TH1D(Form("%s, fccd = %.2f mm, dlf = %.2f", det_names[ch].c_str(), fccd, dlf),
                    Form("%s, fccd = 0.65:0.025:2.4 mm, dlf = %.2f", det_names[ch].c_str(), dlf), 200, 0, 200);

            for (int i = 1; i <= h->GetNbinsX(); ++i) {
                h->SetBinContent(i, gerda::ar39_pdf(ch, h->GetBinCenter(i), fccd, dlf));
            }
            h->Draw("l same");
            histos.push_back(h);
        }

        c->Update();
        auto l = new TLine(40, 0, 40, gPad->GetUymax());
        l->SetLineStyle(2);
        l->Draw();

        c->SaveAs(Form("plots/fccd/ch%i-dlf%.2f.pdf", ch, dlf));
        delete c;
        for (auto h : histos) delete h;
        histos.clear();
    }
}

void dlf_sweep(double fccd=1.) {

    for (auto ch : channels) {
        auto c = new TCanvas(Form("ch%i", ch), "c");
        c->cd();
        std::vector<TH1D*> histos;
        // on grid
        for (double dlf = 0.1; dlf <= 1; dlf += 0.1) {

            auto h = new TH1D(Form("%s, fccd = %.2f mm, dlf = %.2f", det_names[ch].c_str(), fccd, dlf),
                    Form("%s, fccd = %.2f mm, dlf = 0.1:0.05:1", det_names[ch].c_str(), fccd), 200, 0, 200);

            for (int i = 1; i <= h->GetNbinsX(); ++i) {
                h->SetBinContent(i, gerda::ar39_pdf(ch, h->GetBinCenter(i), fccd, dlf));
            }
            h->SetLineColor(kRed);
            h->Draw("l same");
            histos.push_back(h);
        }
        // out grid
        for (double dlf = 0.15; dlf <= 1; dlf += 0.1) {

            auto h = new TH1D(Form("%s, fccd = %.2f mm, dlf = %.2f", det_names[ch].c_str(), fccd, dlf),
                    Form("%s, fccd = %.2f mm, dlf = 0.1:0.05:1", det_names[ch].c_str(), fccd), 200, 0, 200);

            for (int i = 1; i <= h->GetNbinsX(); ++i) {
                h->SetBinContent(i, gerda::ar39_pdf(ch, h->GetBinCenter(i), fccd, dlf));
            }
            h->Draw("l same");
            histos.push_back(h);
        }

        c->Update();
        auto l = new TLine(40, 0, 40, gPad->GetUymax());
        l->SetLineStyle(2);
        l->Draw();

        c->SaveAs(Form("plots/dlf/ch%i-fccd%.2fmm.pdf", ch, fccd));
        delete c;
        for (auto h : histos) delete h;
        histos.clear();
    }
}

void fccd_sweep_fine(double dlf=0.2) {

    for (auto ch : channels) {
        auto c = new TCanvas(Form("ch%i", ch), "c");
        c->cd();
        std::vector<TH1D*> histos;

        // out grid
        freec = 10000;
        for (double fccd = 1.005; fccd <= 1.20; fccd += 0.005) {

            if (std::fabs(fccd-1.05) < 1e-5) continue;
            if (std::fabs(fccd-1.10) < 1e-5) continue;
            if (std::fabs(fccd-1.15) < 1e-5) continue;
            if (std::fabs(fccd-1.20) < 1e-5) continue;

            auto h = new TH1D(Form("%s, fccd = %.3f mm, dlf = %.2f", det_names[ch].c_str(), fccd, dlf),
                    Form("%s, fccd = 1.0:0.005:1.2 mm, dlf = %.2f", det_names[ch].c_str(), dlf), 200, 0, 200);

            for (int i = 1; i <= h->GetNbinsX(); ++i) {
                h->SetBinContent(i, gerda::ar39_pdf(ch, h->GetBinCenter(i), fccd, dlf));
            }
            h->GetXaxis()->SetRangeUser(60, 120);
            h->SetLineColor(freec++);
            if (freec == 10009) freec = 10000;
            h->Draw("l same");
            histos.push_back(h);
        }

        // in grid
        for (double fccd = 1.0; fccd <= 1.20; fccd += 0.05) {

            auto h = new TH1D(Form("%s, fccd = %.3f mm, dlf = %.2f", det_names[ch].c_str(), fccd, dlf),
                    Form("%s, fccd = 1.0:0.005:1.2 mm, dlf = %.2f", det_names[ch].c_str(), dlf), 200, 0, 200);

            for (int i = 1; i <= h->GetNbinsX(); ++i) {
                h->SetBinContent(i, gerda::ar39_pdf(ch, h->GetBinCenter(i), fccd, dlf));
            }
            h->SetLineColor(kRed);
            h->Draw("l same");
            histos.push_back(h);
        }

        c->Update();
        auto l = new TLine(40, 0, 40, gPad->GetUymax());
        l->SetLineStyle(2);
        l->Draw();

        c->SaveAs(Form("plots/fccd/ch%i-dlf%.2f-fine.pdf", ch, dlf));
        delete c;
        for (auto h : histos) delete h;
        histos.clear();
    }
}

void test() {

    float f = 20;
    for (int i = 0; i < 9; ++i) new TColor(freec+i, (0+i*f)/255., (51+i*f)/255., (102+i*f)/255., "blue");

    gROOT->SetBatch();

    fccd_sweep(0.2);
    fccd_sweep(0.5);
    fccd_sweep(0.8);

    dlf_sweep(0.7);
    dlf_sweep(1.);
    dlf_sweep(1.5);
    dlf_sweep(2.4);

    fccd_sweep_fine(0.5);
}

void speedtest() {
    // gerda::ar39_pdf(0, 134.15, 1.57, 0.22);
}
