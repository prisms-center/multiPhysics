// cp_unit_test.inc
// Single-material-point test of crystalPlasticity<dim>::calculatePlasticity against MATLAB.
//
// Include INSIDE crystalPlasticity<dim>::getElementalValues(), just before the
// quadrature-point loop:      #include "cp_unit_test.inc"
// Runs only when env var CP_UNIT_TEST is set, runs once, writes cpp_out.txt, exits.
// Run on ONE MPI rank.  Keep the .inc extension so the src/ glob does not compile it alone.
// Needs <map>, <fstream>, <sstream>, <iomanip>, <cstdlib> (add at the top of the file if missing).

if (std::getenv("CP_UNIT_TEST"))
{
  static bool cp_unit_test_done = false;
  if (!cp_unit_test_done)
  {
    cp_unit_test_done = true;
    const unsigned int c = cellID, q = 0;

    // ---------- read inputs written by export_for_cpp.m (row-major) ----------
    std::map<std::string, std::vector<double>> in;
    {
      std::ifstream fin("cp_unit_test_input.txt");
      AssertThrow(fin.good(), dealii::ExcMessage("cp_unit_test_input.txt not found"));
      std::string line;
      while (std::getline(fin, line)) {
        std::istringstream ss(line); std::string key; double v;
        ss >> key; while (ss >> v) in[key].push_back(v);
      }
    }
    auto M3 = [&](const std::string &k) {
      dealii::FullMatrix<double> A(3, 3);
      for (unsigned int i = 0; i < 3; ++i) for (unsigned int j = 0; j < 3; ++j)
        A(i, j) = in.at(k)[3 * i + j];
      return A;
    };
    auto put = [](std::ostream &o, const std::string &k, const dealii::FullMatrix<double> &A) {
      o << k; for (unsigned int i = 0; i < A.m(); ++i) for (unsigned int j = 0; j < A.n(); ++j)
        o << ' ' << std::setprecision(17) << A(i, j);
      o << '\n';
    };
    auto putv = [](std::ostream &o, const std::string &k, const dealii::Vector<double> &v) {
      o << k; for (unsigned int a = 0; a < v.size(); ++a) o << ' ' << std::setprecision(17) << v[a];
      o << '\n';
    };

    // ---------- overwrite the stored state at (c,q) ----------
    // Every array is checked before it is touched: arrays the model never allocated are
    // reported and skipped instead of segfaulting.
    auto ok = [&](const auto &arr, const char *name) {
      const bool good = c < arr.size() && q < arr[c].size();
      std::cout << "  " << name << (good ? ": set" : ": NOT ALLOCATED, skipped") << std::endl;
      return good;
    };
    std::cout << "CP unit test: setting state at cell " << c << ", q " << q << std::endl;

    if (ok(this->Fp_conv, "Fp_conv")) this->Fp_conv[c][q] = M3("FP_t");
    if (ok(this->Fe_conv, "Fe_conv")) this->Fe_conv[c][q] = M3("FE_t");
    if (ok(this->Fp_iter, "Fp_iter")) this->Fp_iter[c][q] = M3("FP_t");
    if (ok(this->Fe_iter, "Fe_iter")) this->Fe_iter[c][q] = M3("FE_t");

    // C++ stores slip + twin systems together (12 + 1 here); MATLAB has only the 12 slip systems.
    // Copy the slip values and make every extra (twin) system effectively inactive.
    const std::vector<double> &s0 = in.at("s_alpha_t");
    auto set_s = [&](auto &arr, const char *name) {
      if (!ok(arr, name)) return;
      const unsigned int n_cpp = arr[c][q].size();
      AssertThrow(s0.size() <= n_cpp,
                  dealii::ExcMessage("MATLAB has more slip systems than the C++ model"));
      for (unsigned int a = 0; a < n_cpp; ++a) arr[c][q][a] = (a < s0.size()) ? s0[a] : 1.0e12;
      std::cout << "    " << s0.size() << " MATLAB systems, " << n_cpp
                << " C++ systems (extra ones set inactive)" << std::endl;
    };
    set_s(this->s_alpha_conv, "s_alpha_conv");
    set_s(this->s_alpha_iter, "s_alpha_iter");

    // Orientation as a Rodrigues vector (same convention as MATLAB rod2mat / odfpoint).
    // calculatePlasticity reads rot_conv (line 232) and rot_original as the parent (line 235).
    const std::vector<double> &r = in.at("rod");
    auto set_rot = [&](auto &arr, const char *name) {
      if (!ok(arr, name)) return;
      if (arr[c][q].size() < 3) { std::cout << "    size " << arr[c][q].size() << ", skipped" << std::endl; return; }
      for (unsigned int i = 0; i < 3; ++i) arr[c][q][i] = r[i];
    };
    set_rot(this->rot_conv,     "rot_conv");
    set_rot(this->rot_iter,     "rot_iter");
    set_rot(this->rotnew_conv,  "rotnew_conv");
    set_rot(this->rotnew_iter,  "rotnew_iter");
    set_rot(this->rot_original, "rot_original");

    // Rate-dependent state (backstresses, signed slip, Ep, Ep_eff_cum): MATLAB has none.
    auto zero_vec = [&](auto &arr, const char *name) {
      if (!ok(arr, name)) return;
      for (unsigned int k = 0; k < arr[c][q].size(); ++k) arr[c][q][k] = 0;
    };
    zero_vec(this->stateVar_conv,     "stateVar_conv");
    zero_vec(this->stateVar_iter,     "stateVar_iter");
    // MATLAB encodes the twin in Rmat and FP_t, so the C++ twin state must be "untwinned".
    zero_vec(this->twinfraction_conv, "twinfraction_conv");
    zero_vec(this->twinfraction_iter, "twinfraction_iter");
    zero_vec(this->TwinFlag_conv,     "TwinFlag_conv");
    zero_vec(this->TwinFlag_iter,     "TwinFlag_iter");
    if (ok(this->twin_conv, "twin_conv")) this->twin_conv[c][q] = 0;
    if (ok(this->twin_iter, "twin_iter")) this->twin_iter[c][q] = 0;

    // ---------- evaluate ----------
    const dealii::FullMatrix<double> F0 = M3("F_tau");
    auto cauchy_to_P = [](const dealii::FullMatrix<double> &F, const dealii::FullMatrix<double> &T) {
      dealii::FullMatrix<double> Finv(3, 3), P(3, 3);
      Finv.invert(F);
      T.mTmult(P, Finv);                 // T * F^{-T}
      P *= F.determinant();              // J T F^{-T}
      return P;
    };

    this->F = F0;
    std::cout << "CP unit test: calling calculatePlasticity" << std::endl;
    this->calculatePlasticity(c, q, 1);
    std::cout << "CP unit test: calculatePlasticity returned" << std::endl;
    const dealii::FullMatrix<double> T0 = this->T;
    const dealii::FullMatrix<double> P0 = cauchy_to_P(F0, T0);

    std::ofstream out("cpp_out.txt");
    put(out,  "P_tau",       P0);         // from Cauchy T (always consistent)
    put(out,  "P_member",    this->P);    // whatever the code itself stored in P
    put(out,  "T_tau",       T0);
    put(out,  "FP_tau",      this->Fp_iter[c][q]);
    put(out,  "FE_tau",      this->Fe_iter[c][q]);
    putv(out, "s_alpha_tau", this->s_alpha_iter[c][q]);
    // updated orientation: write whichever copies exist
    if (c < this->rot_iter.size()    && q < this->rot_iter[c].size())    putv(out, "rot_iter",    this->rot_iter[c][q]);
    if (c < this->rotnew_iter.size() && q < this->rotnew_iter[c].size()) putv(out, "rotnew_iter", this->rotnew_iter[c][q]);

    // Analytical tangent from the code: Tensor<4,3> dP_dF[i][j][k][l] = dP_ij/dF_kl
    // -> row-major 9x9, same layout as fdperturb_cp3D.
    dealii::FullMatrix<double> dP_an(9, 9);
    for (unsigned int i = 0; i < 3; ++i) for (unsigned int j = 0; j < 3; ++j)
      for (unsigned int k = 0; k < 3; ++k) for (unsigned int l = 0; l < 3; ++l)
        dP_an(3 * i + j, 3 * k + l) = this->dP_dF[i][j][k][l];
    put(out, "dPdF_an", dP_an);

    std::cout << "CP unit test: computing finite-difference tangent" << std::endl;
    // Forward-difference tangent, identical to fdperturb_cp3D (delta = 1e-6, row-major).
    const double delta = 1e-6;
    dealii::FullMatrix<double> dP_fd(9, 9);
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j) {
        dealii::FullMatrix<double> Fp = F0; Fp(i, j) += delta;
        this->F = Fp;
        this->calculatePlasticity(c, q, 0);   // converged state is untouched -> repeatable
        const dealii::FullMatrix<double> Pp = cauchy_to_P(Fp, this->T);
        for (unsigned int k = 0; k < 3; ++k) for (unsigned int l = 0; l < 3; ++l)
          dP_fd(3 * k + l, 3 * i + j) = (Pp(k, l) - P0(k, l)) / delta;
      }
    put(out, "dPdF_fd", dP_fd);

    out.close();
    std::cout << "CP unit test done -> cpp_out.txt" << std::endl;
    std::exit(0);
  }
}
