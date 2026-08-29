/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.Lemma37A
import ABC3.Found.GaloisRep.NorthcottHtJ
import ABC3.Found.GenEll.Elementary
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★`Lemma 3.7` の第 3 の主張・条件 (a) の側（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.18。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

## ★★★★★★★★★★★★★★★★★★これは何か

原文 p.18–19 の**条件 (a) の側の例外集合**の議論を、そのまま形式化する:

> Now suppose that condition (a) is satisfied. Let `C₀` be as in Lemma 3.6. …
> by Lemma 3.6 …, we conclude that `l ≥ 56d·log(l)/log(2)`, hence that
> `(l/14)·deg∞ ≤ ht^Falt + l·log(2)/(28d) + C′` (‡) … Now since `E_L` has at least one
> prime of bad reduction, it follows that `log(2) ≤ d·deg∞`. … we obtain that
> `l·log(2)/(28d) ≤ ht^Falt + C′`. On the other hand, `log(2)/28 ≥ 2/100`, and, by
> assumption, `l ≥ 100d·ht^Falt`, so … **`ht^Falt([E_L]) ≤ C′`**. But this implies,
> by Proposition 3.4, that `[E_L]` belongs to some [fixed] **finite** exceptional set `Exc_d`

★★**その最後の「finite exceptional set」まで取る**——`Proposition 3.4` の有限性は
`§9-1005`（`finite_j_of_htFalt_le`、**無条件**）が与える。

## ★入力

☆本ファイルが受けるのは **`(†)`**（`Lemma 3.5` の結論）だけである。

    `(l/14)·deg∞(E) ≤ ht^Falt(E) + 2·log(l) + C′`

★`Lemma 3.5` に残る唯一の入力は Faltings 高さの同種写像評価であり、
`Found/GaloisRep/Lemma35Concrete.lean` に型で固定してある。
★★**それが landing すれば、本ファイルの `hdag` はそのまま埋まる。**

## ★★他の材料はすべて手元にある

| 材料 | 出どころ |
|---|---|
| `log 2 ≤ d·deg∞`（悪い素点が 1 つある） | ★`minDeltaExp_log_two_le`（`§9-1010`） |
| `Lemma 3.6`（初等的な評価） | ★`Found/GenEll/Elementary.lean`（§1 で計上済み） |
| `ht^Falt` は下に有界 | ★`exists_htFalt_bddBelow`（`§9-1010`） |
| `Proposition 3.4` の有限性 | ★`finite_j_of_htFalt_le`（`§9-1005`、**無条件**） |

## ★定数の取り方

原文の「if we choose `C` sufficiently large [cf. Proposition 3.4] that the inequality of
condition (a) implies that `l ≥ C₀(56d/log 2)^{1+ϵ}`」を、

    `C ≔ |B| + C₀·(56/0.69)^{1+ϵ} + 1`

として明示した（`B` は `ht^Falt` の下界）。★`condA_gives_lemma36` がその計算である。
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve ABC3.Found.GenEll
open scoped Classical

/-! ## ★★★★★★★★数値の核（原文 p.18–19） -/

/-- ★★★★★★★★**原文 p.18–19 の (†) → (‡) → `ht^Falt ≤ C′` の数値の核**。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

★`Ll` は `log(l)`、`L2` は `log(2)`。★★実数だけの主張なので曲線の型は要らない。 -/
theorem lemma_3_7_c_numeric (d D F l Ll L2 C' : ℝ)
    (hd1 : 1 ≤ d) (_hDnn : 0 ≤ D) (hl1 : 1 ≤ l) (_hLlnn : 0 ≤ Ll)
    (hL2lo : 0.69 ≤ L2) (_hL2hi : L2 ≤ 0.70) (hC'0 : 0 ≤ C')
    (hbad : L2 ≤ d * D)
    (hdag : (l / 14) * D ≤ F + 2 * Ll + C')
    (h36 : (56 * d / L2) * Ll ≤ l)
    (ha : 100 * d * F ≤ l) :
    F ≤ C' := by
  have hd0 : (0:ℝ) < d := by linarith
  have hL20 : (0:ℝ) < L2 := by linarith
  have h1 : 2 * Ll * (28 * d) ≤ l * L2 := by
    have h : 56 * d * Ll / L2 ≤ l := by
      rw [div_mul_eq_mul_div] at h36; linarith [h36]
    rw [div_le_iff₀ hL20] at h
    nlinarith [h]
  have h2 : l * L2 ≤ ((l / 14) * D) * (14 * d) := by
    have hl0 : (0:ℝ) ≤ l := by linarith
    nlinarith [mul_le_mul_of_nonneg_left hbad hl0]
  have h3 : l * L2 ≤ (F + C') * (28 * d) := by
    nlinarith [h1, h2, hdag]
  rcases le_or_gt F 0 with hF | hF
  · linarith
  · have hkey : 2 * F * (28 * d) ≤ l * L2 := by nlinarith [ha, hL2lo, hF, hd0]
    nlinarith [h3, hkey, hd0]

/-! ## ★★★★★`C` の取り方——条件 (a) が `Lemma 3.6` の仮説を満たすこと -/

/-- ★★★★★**`C ≔ |B| + C₀·(56/0.69)^{1+ϵ} + 1` なら条件 (a) は `Lemma 3.6` の仮説を含む**。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

★原文「if we choose `C` sufficiently large … that the inequality of condition (a) implies
that `l ≥ C₀(56d/log 2)^{1+ϵ}`」の中身である。
★★`d ≥ 1` から `d·d^ϵ = d^{1+ϵ} ≥ d`、`L2 ≥ 0.69` から `56/L2 ≤ 56/0.69` を使う。 -/
theorem condA_gives_lemma36 (d eps C₀ L2 B F : ℝ) (hd1 : 1 ≤ d) (heps : 0 < eps) (hC₀ : 0 < C₀)
    (hL2lo : 0.69 ≤ L2) (_hL2hi : L2 ≤ 0.70) (hBF : B ≤ F) :
    C₀ * (56 * d / L2) ^ (1 + eps)
      ≤ 100 * d * (F + (|B| + C₀ * ((56:ℝ)/0.69) ^ (1 + eps) + 1) * d ^ eps) := by
  have hd0 : (0:ℝ) < d := by linarith
  have hL20 : (0:ℝ) < L2 := by linarith
  have h1e : (0:ℝ) < 1 + eps := by linarith
  set P : ℝ := d ^ (1 + eps) with hP
  set R : ℝ := ((56:ℝ)/0.69) ^ (1 + eps) with hR
  have hdd : d * d ^ eps = P := by rw [hP, Real.rpow_add hd0, Real.rpow_one]
  have hdle : d ≤ P := by
    rw [hP]
    calc d = d ^ (1:ℝ) := (Real.rpow_one d).symm
      _ ≤ d ^ (1 + eps) := Real.rpow_le_rpow_of_exponent_le hd1 (by linarith)
  have hP0 : (0:ℝ) < P := by rw [hP]; exact Real.rpow_pos_of_pos hd0 _
  have hsplit : (56 * d / L2) ^ (1 + eps) = ((56:ℝ)/L2) ^ (1 + eps) * P := by
    rw [hP, show (56 * d / L2 : ℝ) = (56 / L2) * d by ring]
    exact Real.mul_rpow (by positivity) hd0.le
  have hratio : ((56:ℝ)/L2) ^ (1 + eps) ≤ R := by
    rw [hR]
    exact Real.rpow_le_rpow (by positivity)
      (div_le_div_of_nonneg_left (by norm_num) (by norm_num) hL2lo) h1e.le
  have hR0 : (0:ℝ) < R := by rw [hR]; exact Real.rpow_pos_of_pos (by norm_num) _
  have hBnn : (0:ℝ) ≤ |B| := abs_nonneg B
  have hBge : -|B| ≤ B := neg_abs_le B
  have hexp : 100 * d * (F + (|B| + C₀ * R + 1) * d ^ eps)
      = 100 * d * F + 100 * (|B| + C₀ * R + 1) * P := by rw [← hdd]; ring
  have hlhs : C₀ * (((56:ℝ)/L2) ^ (1 + eps) * P) ≤ C₀ * (R * P) := by
    have h := mul_le_mul_of_nonneg_right hratio hP0.le
    nlinarith [h, hC₀]
  have hdF : -(100 * |B| * P) ≤ 100 * d * F := by
    have h1 : 100 * d * B ≤ 100 * d * F := by nlinarith
    have h2 : -(100 * d * |B|) ≤ 100 * d * B := by nlinarith
    have h3 : -(100 * |B| * P) ≤ -(100 * d * |B|) := by nlinarith
    linarith
  have hE1 : 100 * (|B| + C₀ * R + 1) * P
      = 100 * (|B| * P) + 100 * (C₀ * (R * P)) + 100 * P := by ring
  have hCRP : (0:ℝ) ≤ C₀ * (R * P) := mul_nonneg hC₀.le (mul_nonneg hR0.le hP0.le)
  rw [hsplit, hexp, hE1]
  linarith [hlhs, hdF, hP0, hCRP]

/-! ## ★★★★★★★★★★★★曲線の水準——`ht^Falt ≤ C′` -/

/-- ★★★★★★★★★★★★**条件 (a) ＋ (†) なら `ht^Falt ≤ C′`**。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

☆受けているのは `hdag`（＝ `Lemma 3.5` の結論 (†)）だけである。
★★他の材料——`log 2 ≤ d·deg∞`・`Lemma 3.6`・`ht^Falt` の下界——は**すべて手元にある**。 -/
theorem htFalt_le_of_condA (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, 0 < C ∧ ∀ (C' : ℝ), 0 ≤ C' →
      ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L) [E.IsElliptic]
        (l : ℕ) (p : HeightOneSpectrum (𝓞 L)),
        1 ≤ l →
        minDeltaExp p E ≠ 0 →
        100 * (Module.finrank ℚ L : ℝ)
            * (htFaltOf L E + C * (Module.finrank ℚ L : ℝ) ^ eps) ≤ (l : ℝ) →
        ((l : ℝ) / 14) * degInfOf L E ≤ htFaltOf L E + 2 * Real.log l + C' →
        htFaltOf L E ≤ C' := by
  obtain ⟨C₀, hC₀0, hC₀⟩ := ABC3.Found.GenEll.lemma_3_6 eps heps
  obtain ⟨B, hB⟩ := exists_htFalt_bddBelow
  have hL2lo : (0.69 : ℝ) ≤ Real.log 2 := by linarith [Real.log_two_gt_d9]
  have hL2hi : Real.log 2 ≤ (0.70 : ℝ) := by linarith [Real.log_two_lt_d9]
  refine ⟨|B| + C₀ * ((56:ℝ)/0.69) ^ (1 + eps) + 1, by positivity, ?_⟩
  intro C' hC'0 L _ _ E _ l p hl1 hne hcondA hdag
  set d : ℝ := (Module.finrank ℚ L : ℝ) with hd
  have hd1 : (1:ℝ) ≤ d := by rw [hd]; exact_mod_cast Module.finrank_pos
  have hlR : (1:ℝ) ≤ (l : ℝ) := by exact_mod_cast hl1
  have hLlnn : (0:ℝ) ≤ Real.log l := Real.log_nonneg hlR
  have hΔ : E.Δ ≠ 0 := E.isUnit_Δ.ne_zero
  have hbad : Real.log 2 ≤ d * degInfOf L E := by
    have h := minDeltaExp_log_two_le E hΔ p
    have hge1 : (1:ℤ) ≤ minDeltaExp p E := by
      have := minDeltaExp_nonneg p E; omega
    have hge1R : (1:ℝ) ≤ (minDeltaExp p E : ℝ) := by exact_mod_cast hge1
    nlinarith [h, hge1R, hL2lo]
  have hy1 : (1:ℝ) ≤ 56 * d / Real.log 2 := by
    rw [le_div_iff₀ (by linarith)]
    nlinarith
  have hxge : C₀ * (56 * d / Real.log 2) ^ (1 + eps) ≤ (l : ℝ) :=
    le_trans (condA_gives_lemma36 d eps C₀ (Real.log 2) B (htFaltOf L E) hd1 heps hC₀0
      hL2lo hL2hi (hB L E)) hcondA
  have h36 : (56 * d / Real.log 2) * Real.log l ≤ (l : ℝ) := hC₀ _ _ hy1 hxge
  have hP1 : (1:ℝ) ≤ d ^ eps := Real.one_le_rpow hd1 heps.le
  have hKpos : (0:ℝ) < |B| + C₀ * ((56:ℝ)/0.69) ^ (1 + eps) + 1 := by positivity
  have hdeps : (0:ℝ) < d ^ eps := lt_of_lt_of_le zero_lt_one hP1
  have ha : 100 * d * htFaltOf L E ≤ (l : ℝ) := by
    refine le_trans ?_ hcondA
    have hexp : 100 * d * (htFaltOf L E
          + (|B| + C₀ * ((56:ℝ)/0.69) ^ (1 + eps) + 1) * d ^ eps)
        = 100 * d * htFaltOf L E
          + 100 * d * ((|B| + C₀ * ((56:ℝ)/0.69) ^ (1 + eps) + 1) * d ^ eps) := by ring
    rw [hexp]
    have hnn : (0:ℝ) ≤ 100 * d * ((|B| + C₀ * ((56:ℝ)/0.69) ^ (1 + eps) + 1) * d ^ eps) :=
      mul_nonneg (by linarith) (mul_nonneg hKpos.le hdeps.le)
    linarith
  exact lemma_3_7_c_numeric d (degInfOf L E) (htFaltOf L E) l (Real.log l) (Real.log 2) C'
    hd1 (degInfOf_nonneg E) hlR hLlnn hL2lo hL2hi hC'0 hbad hdag h36 ha

/-! ## ★★★★★★★★★★★★★★★★★★例外集合が有限であること -/

/-- ★★★★★★★★★★★★★★★★★★**条件 (a) の側の例外集合は有限**。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

原文の「But this implies, by Proposition 3.4, that `[E_L]` belongs to some [fixed]
**finite** exceptional set `Exc_d`」までを取る。

★★★次数 `d` を止めた族について有限なので、`d` を動かせば
`Example 1.3, (i)` の意味で **Galois-finite** である。

☆受けているのは `hdag`（＝ `Lemma 3.5` の結論）だけである。 -/
theorem finite_j_of_condA (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, 0 < C ∧ ∀ (C' : ℝ), 0 ≤ C' → ∀ (d : ℕ) {P : Type}
      (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
      (hdeg : ∀ p, Module.finrank ℚ (fld p) ≤ d)
      (E : ∀ p, WeierstrassCurve (fld p)) (hE : ∀ p, (E p).IsElliptic)
      (l : P → ℕ) (q : ∀ p, haveI := hnf p; HeightOneSpectrum (𝓞 (fld p))),
      (∀ p, 1 ≤ l p) →
      (∀ p, haveI := hnf p; haveI := hE p; minDeltaExp (q p) (E p) ≠ 0) →
      (∀ p, haveI := hnf p; haveI := hE p;
        100 * (Module.finrank ℚ (fld p) : ℝ)
          * (htFaltOf (fld p) (E p) + C * (Module.finrank ℚ (fld p) : ℝ) ^ eps) ≤ (l p : ℝ)) →
      (∀ p, haveI := hnf p; haveI := hE p;
        ((l p : ℝ) / 14) * degInfOf (fld p) (E p)
          ≤ htFaltOf (fld p) (E p) + 2 * Real.log (l p) + C') →
      ((fun p : P => haveI := hE p; (((E p).j : fld p) : ℂ)) '' Set.univ).Finite := by
  obtain ⟨C, hCpos, hC⟩ := htFalt_le_of_condA eps heps
  refine ⟨C, hCpos, fun C' hC'0 d P fld hnf hdeg E hE l q hl1 hne hcondA hdag => ?_⟩
  refine Set.Finite.subset
    (finite_j_of_htFalt_le d fld hnf hdeg E hE
      (fun p => haveI := hnf p; haveI := hE p; htFaltOf (fld p) (E p)) (fun p => rfl) C') ?_
  refine Set.image_mono (fun p _ => ?_)
  haveI := hnf p; haveI := hE p
  exact hC C' hC'0 (fld p) (E p) (l p) (q p) (hl1 p) (hne p) (hcondA p) (hdag p)

/-! ## ★出典の紐付け(`.src`)——★★**条つき（`hdag` を受けている）** -/

def lemma_3_7_c_numeric.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(第 3 の主張・条件 (a)——数値の核 (†)→(‡)→ht^Falt ≤ C′)",
    sectionId := "genell-lemma-3-7" }

def condA_gives_lemma36.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(C の取り方——条件 (a) が Lemma 3.6 の仮説を含むこと)",
    sectionId := "genell-lemma-3-7" }

def finite_j_of_condA.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(第 3 の主張・条件 (a) の側の例外集合は有限。★(†) を受ける)",
    sectionId := "genell-lemma-3-7" }

def finite_j_of_condA.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_6(初等的な評価——§1 で計上済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.lemma_3_6") 2,
    .citation "[ABC3]" "minDeltaExp_log_two_le(log 2 ≤ d·deg∞、§9-1010)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.minDeltaExp_log_two_le") 3,
    .citation "[ABC3]" "exists_htFalt_bddBelow(ht^Falt は下に有界、§9-1010)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_htFalt_bddBelow") 2,
    .citation "[ABC3]" "finite_j_of_htFalt_le(Prop 3.4 の有限性——無条件、§9-1005)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.finite_j_of_htFalt_le") 4,
    .otherPaper "[GenEll]"
      ("Lemma 3.5 の結論 (†)——☆本ファイルが受けている**唯一の入力** `hdag` である。" ++
       "★Lemma 3.5 に残る入力は Faltings 高さの同種写像評価だけであり、" ++
       "Found/GaloisRep/Lemma35Concrete.lean に型で固定してある。" ++
       "★★それが landing すれば hdag はそのまま埋まる") 9,
    .implicitStep
      ("★★★★★★★★★★到達点(2026-08-29、第 569): 原文 p.18-19 の**条件 (a) の側の" ++
       "例外集合の議論を最後(「finite exceptional set」)まで取った**。" ++
       "★受けているのは (†) だけで、他の材料——log 2 ≤ d·deg∞(§9-1010)・" ++
       "Lemma 3.6(§1 計上済み)・ht^Falt の下界(§9-1010)・" ++
       "Prop 3.4 の有限性(§9-1005、無条件)——は**すべて手元にある**。" ++
       "★★定数は C = |B| + C₀·(56/0.69)^{1+ϵ} + 1 と明示した" ++
       "(原文の「if we choose C sufficiently large」)") 9,
    .folklore
      ("☆条件 (b) の側の例外集合は Found/GaloisRep/NoMultRedExc.lean(§9-1014)で" ++
       "**無条件に**取れている。★したがって Lemma 3.7 に残るのは" ++
       "**(†) すなわち Lemma 3.5 だけ**である") 8 ]

end ABC3.Found.GaloisRep
