/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInfBaseChange
import ABC3.Found.GaloisRep.DegInfTateParam
import ABC3.Found.GenEll.JScale
import ABC3.Found.GaloisRep.HtFaltJ
import ABC3.Found.GaloisRep.TateParamMap
import ABC3.Found.GaloisRep.TateSetupDvr

/-!
# 第 954 ブロック —— **★★★★★★★★★★★★★★★★悪い素点の局所データを取り出す**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——(D3) の (b) の半分

`minDeltaExp_eq_mul_of_veluMu`（第 904）は `C`・`hC`・`hc4ne`・`hc4`、
すなわち**極小モデルとその `c₄` が単元であること**を受ける。

★`SemistableAt p W` の定義は

    `minDeltaExp p W = 0` ∨ `∃ C, IsMinimal (primeSubring p) (C • W) ∧ ∃ h, valAdd p c₄ = 0`

だから、**悪い素点（`jExp p W < 0`）では左の選択肢が落ちる**——
`minDeltaExp p W = max 0 (-jExp p W) = -jExp p W > 0` だからである。

☆したがって半安定性だけから 4 つのデータがそのまま出る。

| 定理 | 内容 |
|---|---|
| `minDeltaExp_pos_of_jExp_neg` | ★`jExp < 0` なら `Δ_min > 0` |
| `exists_minimal_c4_unit_of_jExp_neg` | ★★★★★★★★★★★★★★★★**極小モデルと `c₄` の単元性** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField
open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-- ★**半安定で `jExp < 0` なら `Δ_min > 0`**。 -/
theorem minDeltaExp_pos_of_jExp_neg (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [W.IsElliptic] (hss : SemistableAt p W)
    (hj : jExp p W < 0) : 0 < minDeltaExp p W := by
  rw [minDeltaExp_eq_maxJ_of_semistable p W hss]
  omega

/-- ★★★★★★★★★★★★★★★★**悪い素点では半安定性が
極小モデルと `c₄` の単元性を直に与える**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 954）**——これが
`minDeltaExp_eq_mul_of_veluMu`（第 904）の `C`・`hC`・`hc4ne`・`hc4` の出どころである。
☆`SemistableAt` の左の選択肢（`Δ_min = 0`）は `jExp < 0` と矛盾する。 -/
theorem exists_minimal_c4_unit_of_jExp_neg (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [W.IsElliptic] (hss : SemistableAt p W)
    (hj : jExp p W < 0) :
    ∃ C : WeierstrassCurve.VariableChange L, IsMinimal (primeSubring p) (C • W) ∧
      ∃ h : (C • W).c₄ ≠ 0, valAdd p (Units.mk0 ((C • W).c₄) h) = 0 := by
  rcases hss with h0 | hdata
  · exact absurd h0 (minDeltaExp_pos_of_jExp_neg p W (Or.inl h0) hj).ne'
  · exact hdata

/-! ## ★★★★★★★★★★★★第 973 —— `v_p(Δ) = −jExp`（極小モデルで）

★第 909（`hasMultiplicativeReduction_baseChange`）は `0 < v_p(Δ)` を受ける。
☆極小モデルでは `v_p(c₄) = 0` なので、`j = Δ⁻¹c₄³` より `v_p(j) = −v_p(Δ)`、
すなわち `v_p(Δ) = −jExp`。★悪い素点（`jExp < 0`）ではこれが正である。 -/

/-- ★★★★★★★★★★★★**極小モデルでは `v_p(Δ) = −jExp`**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`j = Δ⁻¹c₄³` を単元の等式に直し、`valAdd` の乗法性で割るだけである。 -/
theorem valAdd_Delta_eq_neg_jExp (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [W.IsElliptic]
    (hΔ : W.Δ ≠ 0) (hc4ne : W.c₄ ≠ 0)
    (hc4 : valAdd p (Units.mk0 W.c₄ hc4ne) = 0) :
    valAdd p (Units.mk0 W.Δ hΔ) = - jExp p W := by
  have hjeq : W.j = W.Δ⁻¹ * W.c₄ ^ 3 := ABC3.Found.GenEll.j_eq_inv_Delta_mul W
  have hj : W.j ≠ 0 := by
    rw [hjeq]
    exact mul_ne_zero (inv_ne_zero hΔ) (pow_ne_zero 3 hc4ne)
  have hunit : Units.mk0 W.j hj
      = (Units.mk0 W.c₄ hc4ne) ^ 3 * (Units.mk0 W.Δ hΔ)⁻¹ := by
    ext
    show W.j = ((Units.mk0 W.c₄ hc4ne) ^ 3 * (Units.mk0 W.Δ hΔ)⁻¹ : Lˣ)
    push_cast
    rw [hjeq]
    show W.Δ⁻¹ * W.c₄ ^ 3 = W.c₄ ^ 3 * (W.Δ)⁻¹
    ring
  have hJ : jExp p W = valAdd p (Units.mk0 W.j hj) := dif_neg hj
  rw [hJ, hunit, valAdd_mul, valAdd_pow, valAdd_inv, hc4]
  omega

/-- ★★★★★★★★**悪い素点では極小モデルの `v_p(Δ)` は正**——第 909 の仮説である。 -/
theorem valAdd_Delta_pos_of_jExp_neg (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [W.IsElliptic]
    (hΔ : W.Δ ≠ 0) (hc4ne : W.c₄ ≠ 0)
    (hc4 : valAdd p (Units.mk0 W.c₄ hc4ne) = 0) (hj : jExp p W < 0) :
    0 < valAdd p (Units.mk0 W.Δ hΔ) := by
  rw [valAdd_Delta_eq_neg_jExp p W hΔ hc4ne hc4]
  omega

def valAdd_Delta_eq_neg_jExp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(極小モデルでは v_p(Δ) = −jExp。★無条件)",
    sectionId := "genell-lemma-3-5" }

def valAdd_Delta_pos_of_jExp_neg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点では極小モデルの v_p(Δ) は正。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★第 976 —— 悪い素点で局所の還元型を出す

★第 972 は `E ⊗ Lv` の**極小性**と**分裂乗法還元**を受ける。
☆第 974 で測ったとおり `IsMinimal` は変数変換で不変ではないので、
渡すのは `C • E`（`SemistableAt` が与える極小モデル）である。

★本ブロックは `SemistableAt` の 4 つのデータ（`C`・`hC`・`hc4ne`・`hc4`、第 954）から

* `(C • E) ⊗ Lv` の極小性（第 908）
* `(C • E) ⊗ Lv` の**乗法**還元（第 909、`0 < v_p(Δ)` は第 973 が与える）

を出す。☆残るのは分裂性だけで、それは第 963 の二者択一が扱う。 -/

section BadPrimeLocal

variable {Lv : Type} [Field Lv] [Algebra L Lv]
  {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [Algebra R Lv] [IsFractionRing R Lv]

/-- ★**`jExp` は変数変換で不変**——`j` が不変だから。 -/
theorem jExp_variableChange (p : HeightOneSpectrum (𝓞 L))
    (E : WeierstrassCurve L) [E.IsElliptic] (C : WeierstrassCurve.VariableChange L) :
    jExp p (C • E) = jExp p E :=
  jExp_congr_j p (C • E) E (WeierstrassCurve.variableChange_j E C)

/-- ★★★★★★★★**極小モデルは完備化でも極小**（第 908 を `C • E` に当てた形）。 -/
theorem isMinimal_baseChange_at_bad_prime (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = (HeightOneSpectrum.valuation L p) x)
    (E : WeierstrassCurve L) [E.IsElliptic]
    (C : WeierstrassCurve.VariableChange L)
    (hC : WeierstrassCurve.IsMinimal (primeSubring p) (C • E))
    (hc4ne : (C • E).c₄ ≠ 0) (hc4 : valAdd p (Units.mk0 ((C • E).c₄) hc4ne) = 0) :
    WeierstrassCurve.IsMinimal R ((C • E).baseChange Lv) := by
  haveI := hC
  exact isMinimal_baseChange_of_c4 p hp (C • E)
    (variableChange_Delta_ne_zero E (E.isUnit_Δ.ne_zero) C) hc4ne hc4

/-- ★★★★★★★★★★★★**悪い素点では完備化で乗法還元をもつ**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 976）**——`0 < v_p(Δ)` は第 973 が `jExp < 0` から与える。
☆残るのは分裂性だけで、それは第 963 の二者択一が扱う。 -/
theorem hasMultiplicativeReduction_at_bad_prime (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = (HeightOneSpectrum.valuation L p) x)
    (E : WeierstrassCurve L) [E.IsElliptic]
    (C : WeierstrassCurve.VariableChange L)
    (hC : WeierstrassCurve.IsMinimal (primeSubring p) (C • E))
    (hc4ne : (C • E).c₄ ≠ 0) (hc4 : valAdd p (Units.mk0 ((C • E).c₄) hc4ne) = 0)
    (hj : jExp p E < 0)
    [WeierstrassCurve.IsMinimal R ((C • E).baseChange Lv)] :
    WeierstrassCurve.HasMultiplicativeReduction R ((C • E).baseChange Lv) := by
  haveI := hC
  have hΔ : (C • E).Δ ≠ 0 := variableChange_Delta_ne_zero E (E.isUnit_Δ.ne_zero) C
  have hjC : jExp p (C • E) < 0 := by rw [jExp_variableChange p E C]; exact hj
  exact hasMultiplicativeReduction_baseChange p hp (C • E) hΔ hc4ne hc4
    (valAdd_Delta_pos_of_jExp_neg p (C • E) hΔ hc4ne hc4 hjC)

end BadPrimeLocal

def jExp_variableChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(jExp は変数変換で不変。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isMinimal_baseChange_at_bad_prime.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(極小モデルは完備化でも極小。★無条件)",
    sectionId := "genell-lemma-3-5" }

def hasMultiplicativeReduction_at_bad_prime.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点では完備化で乗法還元をもつ。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★第 978 —— `hcop` の出どころ

★第 972 の `hcop` は **Tate 母数の付値**についての「`l` と互いに素」である。
☆原文の仮定「`l` is prime to the local heights」は `jExp` の言葉なので、
両者を繋ぐ橋が要る——それが第 932（`jExp_eq_neg_vAdd_of_j_tateCurveAt`）である。

★第 932 の仮説はすべて出る:

* Tate モデルの楕円性・`c₄ ≠ 0`・`j` の一致 → `tateModel_baseChange`（第 944）
* `1/j` の評価が `0` でない → `evalAdic_tateJinvSeries_eq_mul_unit`（`q·単元` だから）
* `E.j ≠ 0`・`E.c₄ ≠ 0` → `jExp p E < 0`（`j = 0` なら `jExp = 0`） -/

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★**Tate 母数の付値は `−jExp`**。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★★★**2026-09-01（第 978）**——第 932 を `q = tateParamR (E ⊗ Lv) h` に当てた形。
☆`mkTateSetup` の `v`・`Q` は定義上 `tateDvrVal`・`Units.mk0 (φ q)` なので `rfl` で繋がる。 -/
theorem vAdd_tateParam_eq_neg_jExp {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = (HeightOneSpectrum.valuation L p) x)
    (E : WeierstrassCurve L) [E.IsElliptic]
    [(E.baseChange Lv).IsElliptic] [WeierstrassCurve.IsMinimal R (E.baseChange Lv)]
    (h : WeierstrassCurve.HasSplitMultiplicativeReduction R (E.baseChange Lv))
    (hj : jExp p E < 0) :
    vAdd (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h)
        (tateParamR_mem (E.baseChange Lv) h) (tateParamR_ne_zero (E.baseChange Lv) h)).v
      (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h)
        (tateParamR_mem (E.baseChange Lv) h)
        (tateParamR_ne_zero (E.baseChange Lv) h)).Q = - jExp p E := by
  have hq := tateParamR_mem (E.baseChange Lv) h
  have hq0 := tateParamR_ne_zero (E.baseChange Lv) h
  obtain ⟨hq', C₀, hne, hCE⟩ := tateParamR_spec (E.baseChange Lv) h
  have hbase := tateModel_baseChange (E.baseChange Lv) h hCE
  haveI : ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map
      (algebraMap R Lv)).IsElliptic := by rw [hbase]; infer_instance
  have hjne : E.j ≠ 0 := by
    intro hc; rw [jExp, dif_pos hc] at hj; omega
  have hc4E : E.c₄ ≠ 0 := by
    intro hc; apply hjne; rw [ABC3.Found.GenEll.j_eq_inv_Delta_mul, hc]; ring
  have hqne : algebraMap R Lv (tateParamR (E.baseChange Lv) h) ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective R Lv)).2 hq0
  obtain ⟨u, hu, hueq⟩ := evalAdic_tateJinvSeries_eq_mul_unit
    (I := IsLocalRing.maximalIdeal R) (tateParamR (E.baseChange Lv) h) hq
  have hev : algebraMap R Lv
      (evalAdic tateJinvSeries (tateParamR (E.baseChange Lv) h) hq) ≠ 0 := by
    rw [hueq, map_mul]
    exact mul_ne_zero hqne ((hu.map (algebraMap R Lv)).ne_zero)
  have hc4Lv : (E.baseChange Lv).c₄ ≠ 0 := by
    show (E.map (algebraMap L Lv)).c₄ ≠ 0
    rw [WeierstrassCurve.map_c₄]
    exact (map_ne_zero_iff _ (algebraMap L Lv).injective).2 hc4E
  have hc4T : algebraMap R Lv
      (tateCurveAt (tateParamR (E.baseChange Lv) h) hq).c₄ ≠ 0 := by
    have hmap : algebraMap R Lv (tateCurveAt (tateParamR (E.baseChange Lv) h) hq).c₄
        = ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map (algebraMap R Lv)).c₄ :=
      (WeierstrassCurve.map_c₄ _ _).symm
    rw [hmap, hbase, WeierstrassCurve.variableChange_c₄]
    exact mul_ne_zero (by simp) hc4Lv
  have hjeq : (E.baseChange Lv).j
      = ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map (algebraMap R Lv)).j := by
    rw [ABC3.Found.GenEll.j_congr_curve hbase, WeierstrassCurve.variableChange_j]
  have hkey := jExp_eq_neg_vAdd_of_j_tateCurveAt p hp E
    (tateParamR (E.baseChange Lv) h) hq hc4T hev hqne hjne hjeq
  rw [hkey, neg_neg]
  rfl

open scoped Classical in
/-- ★★★★★★★★★★★★**原文の「`l` は局所高さと互いに素」を第 972 の `hcop` に直す**。 -/
theorem not_dvd_vAdd_tateParam_of_not_dvd_jExp {Lv : Type} [Field Lv] [CharZero Lv]
    [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = (HeightOneSpectrum.valuation L p) x)
    (E : WeierstrassCurve L) [E.IsElliptic]
    [(E.baseChange Lv).IsElliptic] [WeierstrassCurve.IsMinimal R (E.baseChange Lv)]
    (h : WeierstrassCurve.HasSplitMultiplicativeReduction R (E.baseChange Lv))
    (hj : jExp p E < 0) {l : ℕ} (hcop : ¬ ((l : ℤ) ∣ jExp p E)) :
    ¬ ((l : ℤ) ∣ vAdd (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h)
        (tateParamR_mem (E.baseChange Lv) h) (tateParamR_ne_zero (E.baseChange Lv) h)).v
      (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h)
        (tateParamR_mem (E.baseChange Lv) h)
        (tateParamR_ne_zero (E.baseChange Lv) h)).Q) := by
  rw [vAdd_tateParam_eq_neg_jExp p hp E h hj]
  simpa using hcop

def vAdd_tateParam_eq_neg_jExp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 母数の付値は −jExp。★無条件)",
    sectionId := "genell-def-3-3" }

def not_dvd_vAdd_tateParam_of_not_dvd_jExp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(l は局所高さと互いに素——hcop の形。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★第 995 —— `E′` の側は `j` だけでよい

★★★★**測定**: 第 892／972 は `E′` についても
分裂乗法還元・極小モデル・`C′`・`hc4′` を要求していた。
☆だが結論 `minDeltaExp p E′ = l · minDeltaExp p E` は、半安定性のもとでは

    `minDeltaExp = max(0, −jExp)`（`minDeltaExp_eq_maxJ_of_semistable`、在庫）

なので **`jExp p E′ = l · jExp p E` さえ出れば済む**。

★そして `jExp` は第 932（`jExp_eq_neg_vAdd_of_j_tateCurveAt`）により
**`j` の一致だけ**から出る——`E′` の側に分裂性も極小モデルも要らない。
☆`q^l` の付値は `l · v(q)`（`vAdd_pow`）だから、そのまま `l` 倍になる。 -/

/-- ★★★★★★★★★★★★★★★★★★★★**`E′` の `j` が `E_{q^l}` の `j` なら
`jExp p E′ = l · jExp p E`**。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★★★**2026-09-01（第 995）**——`E′` の側に**分裂性も極小モデルも要らない**。
☆第 932 を `E` と `E′` の両方に当て、`vAdd (q^l) = l · vAdd q` で割るだけである。 -/
theorem jExp_eq_mul_of_tateParam_pow {Lv : Type} [Field Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = (HeightOneSpectrum.valuation L p) x)
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange Lv).IsElliptic] [(E'.baseChange Lv).IsElliptic]
    {l : ℕ} (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R)
    (hql : q ^ l ∈ IsLocalRing.maximalIdeal R)
    [((tateCurveAt q hq).map (algebraMap R Lv)).IsElliptic]
    [((tateCurveAt (q ^ l) hql).map (algebraMap R Lv)).IsElliptic]
    (hc4 : algebraMap R Lv (tateCurveAt q hq).c₄ ≠ 0)
    (hev : algebraMap R Lv (evalAdic tateJinvSeries q hq) ≠ 0)
    (hc4' : algebraMap R Lv (tateCurveAt (q ^ l) hql).c₄ ≠ 0)
    (hev' : algebraMap R Lv (evalAdic tateJinvSeries (q ^ l) hql) ≠ 0)
    (hqne : algebraMap R Lv q ≠ 0) (hqlne : algebraMap R Lv (q ^ l) ≠ 0)
    (hjE : E.j ≠ 0) (hjE' : E'.j ≠ 0)
    (heq : (E.baseChange Lv).j = ((tateCurveAt q hq).map (algebraMap R Lv)).j)
    (heq' : (E'.baseChange Lv).j
      = ((tateCurveAt (q ^ l) hql).map (algebraMap R Lv)).j) :
    jExp p E' = (l : ℤ) * jExp p E := by
  rw [jExp_eq_neg_vAdd_of_j_tateCurveAt p hp E q hq hc4 hev hqne hjE heq,
    jExp_eq_neg_vAdd_of_j_tateCurveAt p hp E' (q ^ l) hql hc4' hev' hqlne hjE' heq']
  have hpow : (Units.mk0 (algebraMap R Lv (q ^ l)) hqlne)
      = (Units.mk0 (algebraMap R Lv q) hqne) ^ l := by
    ext
    simp only [Units.val_mk0, Units.val_pow_eq_pow_val, map_pow]
  rw [hpow, vAdd_pow]
  ring

/-- ★★★★★★★★★★★★**`jExp` が `l` 倍なら `Δ_min` も `l` 倍**（半安定・悪い素点で）。

☆`minDeltaExp = max(0, −jExp)`（在庫の `minDeltaExp_eq_maxJ_of_semistable`）だから、
`jExp p E < 0` なら両辺とも `−jExp` の側が効く。 -/
theorem minDeltaExp_eq_mul_of_jExp_mul (p : HeightOneSpectrum (𝓞 L))
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (hss : SemistableAt p E) (hss' : SemistableAt p E')
    {l : ℕ} (hj : jExp p E < 0) (hmul : jExp p E' = (l : ℤ) * jExp p E) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  rw [minDeltaExp_eq_maxJ_of_semistable p E hss,
    minDeltaExp_eq_maxJ_of_semistable p E' hss', hmul]
  have hl : (0 : ℤ) ≤ (l : ℤ) := Int.natCast_nonneg l
  rcases eq_or_lt_of_le hl with h0 | hpos
  · rw [← h0]; simp
  · rw [max_eq_right (by nlinarith), max_eq_right (by omega)]
    ring

/-! ## ★★★★★★★★★★★★第 1001 —— `SemistableAt` は変数変換で不変

★第 973／976 が与えるのは `E` ではなく **`C • E`** についての
極小性・乗法還元である（`C` は `p` で極小にする変数変換）。
☆したがって第 1000 は `C • E` に当てることになる。

★そのとき `SemistableAt p (C • E)` が要るが、`SemistableAt` の定義は
「`Δ_min = 0`」または「ある `C₀` で極小かつ `c₄` が単元」なので、
**変数変換で不変**である——`C₀` を `C₀ * C⁻¹` に取り替えるだけ。 -/

/-- ★★★★★★★★★★★★**`SemistableAt` は変数変換で不変**。 -/
theorem semistableAt_variableChange (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) (C : WeierstrassCurve.VariableChange L)
    (h : SemistableAt p W) : SemistableAt p (C • W) := by
  rcases h with h0 | ⟨C₀, hC₀, hc4, hv⟩
  · exact Or.inl (by rw [minDeltaExp_variableChange p W C]; exact h0)
  · refine Or.inr ⟨C₀ * C⁻¹, ?_, ?_⟩
    · rw [mul_smul, inv_smul_smul]; exact hC₀
    · rw [mul_smul, inv_smul_smul]; exact ⟨hc4, hv⟩

def semistableAt_variableChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.5(SemistableAt は変数変換で不変。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★第 998 —— `E′.j ≠ 0` は `j` の一致から出る

★第 997 は `hjE′ : E′.j ≠ 0` を仮説で受けていた。
☆だが `j(E_{q₀} ⊗ Lv) = (1/j の評価)⁻¹` であり（`j_tateCurveAt_inv`）、
`1/j` の評価は `q₀ · 単元` だから `q₀ ≠ 0` なら `0` でない。
★したがって **`j` の一致からそのまま出る**——仮説を 1 本減らせる。 -/

/-- ★★★★★★★★★★★★**`j` が Tate 曲線の `j` に等しければ `j ≠ 0`**。 -/
theorem j_ne_zero_of_j_tateCurveAt {Lv : Type} [Field Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (E' : WeierstrassCurve L) [E'.IsElliptic] [(E'.baseChange Lv).IsElliptic]
    (q₀ : R) (hq₀ : q₀ ∈ IsLocalRing.maximalIdeal R) (hq₀0 : q₀ ≠ 0)
    [((tateCurveAt q₀ hq₀).map (algebraMap R Lv)).IsElliptic]
    (hjj : (E'.baseChange Lv).j = ((tateCurveAt q₀ hq₀).map (algebraMap R Lv)).j) :
    E'.j ≠ 0 := by
  have hq₀ne : algebraMap R Lv q₀ ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective R Lv)).2 hq₀0
  have hc4T : algebraMap R Lv (tateCurveAt q₀ hq₀).c₄ ≠ 0 :=
    ((tateCurveAt_c4_isUnit q₀ hq₀).map (algebraMap R Lv)).ne_zero
  obtain ⟨u, hu, hueq⟩ := evalAdic_tateJinvSeries_eq_mul_unit
    (I := IsLocalRing.maximalIdeal R) q₀ hq₀
  have hev : algebraMap R Lv (evalAdic tateJinvSeries q₀ hq₀) ≠ 0 := by
    rw [hueq, map_mul]
    exact mul_ne_zero hq₀ne ((hu.map (algebraMap R Lv)).ne_zero)
  have hjT : ((tateCurveAt q₀ hq₀).map (algebraMap R Lv)).j ≠ 0 := by
    rw [j_tateCurveAt_inv q₀ hq₀ hc4T]
    exact inv_ne_zero hev
  have hjmap : (E'.baseChange Lv).j = algebraMap L Lv E'.j := by
    have hΔbc : (E'.baseChange Lv).Δ = algebraMap L Lv E'.Δ := WeierstrassCurve.map_Δ _ _
    have hcbc : (E'.baseChange Lv).c₄ = algebraMap L Lv E'.c₄ := WeierstrassCurve.map_c₄ _ _
    rw [ABC3.Found.GenEll.j_eq_inv_Delta_mul, ABC3.Found.GenEll.j_eq_inv_Delta_mul,
      hΔbc, hcbc, map_mul, map_inv₀, map_pow]
  intro hc
  apply hjT
  rw [← hjj, hjmap, hc, map_zero]

def j_ne_zero_of_j_tateCurveAt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(j が Tate 曲線の j なら 0 でない。★無条件)",
    sectionId := "genell-def-3-3" }

/-! ## ★★★★★★★★★★★★★★★★★★★★第 997 —— `q` を `tateParamR` に固定した形

★第 995（`jExp_eq_mul_of_tateParam_pow`）は `q`・`hc4`・`hev`・`hqne` を素で受ける。
☆本ブロックは `q := tateParamR (E ⊗ Lv) h` に固定して**それらをすべて自前で作る**。
第 978（`vAdd_tateParam_eq_neg_jExp`）の準備がそのまま使える。

★★結果、`E′` の側に残る仮説は

* `E′.j ≠ 0`
* `j(E′ ⊗ Lv) = j(E_{q^l} ⊗ Lv)`

の 2 本だけになる——**分裂性も極小モデルも `C′` も要らない**。 -/

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★**`E′` の `j` が `E_{q_E^l}` の `j` なら
`jExp p E′ = l · jExp p E`**（`q` を `tateParamR` に固定した形）。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★★★**2026-09-01（第 997）**——第 995 の仮説をすべて `E` の側の
分裂乗法還元から作る。☆第 978 の準備の再利用である。 -/
theorem jExp_eq_mul_of_j_tate_pow {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = (HeightOneSpectrum.valuation L p) x)
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange Lv).IsElliptic] [(E'.baseChange Lv).IsElliptic]
    [WeierstrassCurve.IsMinimal R (E.baseChange Lv)]
    (h : WeierstrassCurve.HasSplitMultiplicativeReduction R (E.baseChange Lv))
    (hj : jExp p E < 0)
    {l : ℕ}
    (hql : (tateParamR (E.baseChange Lv) h) ^ l ∈ IsLocalRing.maximalIdeal R)
    [((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).IsElliptic]
    (hjj : (E'.baseChange Lv).j
      = ((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
        (algebraMap R Lv)).j) :
    jExp p E' = (l : ℤ) * jExp p E := by
  have hq := tateParamR_mem (E.baseChange Lv) h
  have hq0 := tateParamR_ne_zero (E.baseChange Lv) h
  obtain ⟨hq', C₀, hne, hCE⟩ := tateParamR_spec (E.baseChange Lv) h
  have hbase := tateModel_baseChange (E.baseChange Lv) h hCE
  haveI : ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map
      (algebraMap R Lv)).IsElliptic := by rw [hbase]; infer_instance
  -- ★`E` の側（第 978 と同じ）
  have hjne : E.j ≠ 0 := by
    intro hc; rw [jExp, dif_pos hc] at hj; omega
  have hc4E : E.c₄ ≠ 0 := by
    intro hc; apply hjne; rw [ABC3.Found.GenEll.j_eq_inv_Delta_mul, hc]; ring
  have hqne : algebraMap R Lv (tateParamR (E.baseChange Lv) h) ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective R Lv)).2 hq0
  obtain ⟨u, hu, hueq⟩ := evalAdic_tateJinvSeries_eq_mul_unit
    (I := IsLocalRing.maximalIdeal R) (tateParamR (E.baseChange Lv) h) hq
  have hev : algebraMap R Lv
      (evalAdic tateJinvSeries (tateParamR (E.baseChange Lv) h) hq) ≠ 0 := by
    rw [hueq, map_mul]
    exact mul_ne_zero hqne ((hu.map (algebraMap R Lv)).ne_zero)
  have hc4Lv : (E.baseChange Lv).c₄ ≠ 0 := by
    show (E.map (algebraMap L Lv)).c₄ ≠ 0
    rw [WeierstrassCurve.map_c₄]
    exact (map_ne_zero_iff _ (algebraMap L Lv).injective).2 hc4E
  have hc4T : algebraMap R Lv
      (tateCurveAt (tateParamR (E.baseChange Lv) h) hq).c₄ ≠ 0 := by
    have hmap : algebraMap R Lv (tateCurveAt (tateParamR (E.baseChange Lv) h) hq).c₄
        = ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map (algebraMap R Lv)).c₄ :=
      (WeierstrassCurve.map_c₄ _ _).symm
    rw [hmap, hbase, WeierstrassCurve.variableChange_c₄]
    exact mul_ne_zero (by simp) hc4Lv
  have hjeq : (E.baseChange Lv).j
      = ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map (algebraMap R Lv)).j := by
    rw [ABC3.Found.GenEll.j_congr_curve hbase, WeierstrassCurve.variableChange_j]
  -- ☆`q^l` の側——`c₄` は `1` から始まるので常に単元、`1/j` は `q^l · 単元`
  have hqlne : algebraMap R Lv ((tateParamR (E.baseChange Lv) h) ^ l) ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective R Lv)).2 (pow_ne_zero l hq0)
  have hc4T' : algebraMap R Lv
      (tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).c₄ ≠ 0 :=
    ((tateCurveAt_c4_isUnit ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).ne_zero
  obtain ⟨u', hu', hueq'⟩ := evalAdic_tateJinvSeries_eq_mul_unit
    (I := IsLocalRing.maximalIdeal R) ((tateParamR (E.baseChange Lv) h) ^ l) hql
  have hev' : algebraMap R Lv
      (evalAdic tateJinvSeries ((tateParamR (E.baseChange Lv) h) ^ l) hql) ≠ 0 := by
    rw [hueq', map_mul]
    exact mul_ne_zero hqlne ((hu'.map (algebraMap R Lv)).ne_zero)
  have hjE' : E'.j ≠ 0 := j_ne_zero_of_j_tateCurveAt E'
    ((tateParamR (E.baseChange Lv) h) ^ l) hql (pow_ne_zero l hq0) hjj
  exact jExp_eq_mul_of_tateParam_pow p hp E E' (tateParamR (E.baseChange Lv) h) hq hql
    hc4T hev hc4T' hev' hqne hqlne hjne hjE' hjeq hjj

/-- ★★★★★★★★★★★★★★★★**`Δ_min` まで一気に**——半安定性と `j` の一致だけで。 -/
theorem minDeltaExp_eq_mul_of_j_tate_pow {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = (HeightOneSpectrum.valuation L p) x)
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange Lv).IsElliptic] [(E'.baseChange Lv).IsElliptic]
    [WeierstrassCurve.IsMinimal R (E.baseChange Lv)]
    (h : WeierstrassCurve.HasSplitMultiplicativeReduction R (E.baseChange Lv))
    (hss : SemistableAt p E) (hss' : SemistableAt p E')
    (hj : jExp p E < 0)
    {l : ℕ}
    (hql : (tateParamR (E.baseChange Lv) h) ^ l ∈ IsLocalRing.maximalIdeal R)
    [((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).IsElliptic]
    (hjj : (E'.baseChange Lv).j
      = ((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
        (algebraMap R Lv)).j) :
    minDeltaExp p E' = l * minDeltaExp p E :=
  minDeltaExp_eq_mul_of_jExp_mul p E E' hss hss' hj
    (jExp_eq_mul_of_j_tate_pow p hp E E' h hj hql hjj)

def jExp_eq_mul_of_j_tate_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(E′ の j が E_{q_E^l} の j なら jExp は l 倍。★無条件)",
    sectionId := "genell-def-3-3" }

def minDeltaExp_eq_mul_of_j_tate_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(半安定性と j の一致だけで Δ_min が l 倍。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★第 1055 —— `Δ(E_{q^l})` の付値は `l·v(q)`

★第 1054 の節点で要る局所の計算のうち、**Tate 曲線側の一段**である。
☆`Δ(E_q) = q · 単元`（在庫 `tateCurveAt_Delta_eq_mul_unit`）なので、
`q^l` に当てれば `v(Δ(E_{q^l})) = l·v(q)` である。 -/

/-- ★★★★★★★★**`Δ(E_{q^l})` の付値は `l·v(q)`**（第 1055）。 -/
theorem vAdd_Delta_tateCurveAt_pow {R : Type} [CommRing R] [IsDomain R]
    [IsDiscreteValuationRing R] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]
    (q : R) {l : ℕ} (hql : q ^ l ∈ IsLocalRing.maximalIdeal R)
    (hΔ : algebraMap R K (tateCurveAt (q ^ l) hql).Δ ≠ 0)
    (hqne : algebraMap R K q ≠ 0) :
    vAdd (tateDvrVal R K) (Units.mk0 (algebraMap R K (tateCurveAt (q ^ l) hql).Δ) hΔ)
      = l * vAdd (tateDvrVal R K) (Units.mk0 (algebraMap R K q) hqne) := by
  obtain ⟨u, hu, hΔu⟩ := tateCurveAt_Delta_eq_mul_unit (q ^ l) hql
  have hune : algebraMap R K u ≠ 0 := (hu.map (algebraMap R K)).ne_zero
  have hsplit : (Units.mk0 (algebraMap R K (tateCurveAt (q ^ l) hql).Δ) hΔ)
      = (Units.mk0 (algebraMap R K q) hqne) ^ l * (Units.mk0 (algebraMap R K u) hune) := by
    ext
    simp only [Units.val_mk0, Units.val_mul, Units.val_pow_eq_pow_val, hΔu, map_mul, map_pow]
  rw [hsplit, vAdd_mul, vAdd_pow, tateDvrVal_eq_zero_of_isUnit u hu hune, add_zero]

def vAdd_Delta_tateCurveAt_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Δ(E_{q^l}) の付値は l·v(q)。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★出典の紐付け(`.src`) -/

def minDeltaExp_pos_of_jExp_neg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(半安定で jExp < 0 なら Δ_min > 0。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_minimal_c4_unit_of_jExp_neg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点では半安定性が極小モデルと c₄ の単元性を与える。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
