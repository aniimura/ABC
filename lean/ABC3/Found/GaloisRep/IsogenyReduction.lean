/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.HtFaltCovolume
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★同種写像評価は**1 つの局所不等式**に帰着する（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★★★★★★★★★★★★★★★これは何か

§3・§4 の残り 5 項目（`Lemma 3.5` → `Lemma 3.7` → `Theorem 3.8` →
`Corollary 4.3` / `Corollary 4.4`）は、`Lemma 3.5` の唯一の残り入力

    `ht^Falt(E/H) ≤ ht^Falt(E) + 2·log(l)`

を待っている。★★**本ファイルはそれを 1 つの不等式に帰着させる**:

    **`(Σ_σ log‖u′_σ‖ − Σ_σ log‖u_σ‖) − (Σᶠ_p neronExp_p(E′)·log N(p) − Σᶠ_p neronExp_p(E)·log N(p))`**
    **`  ≤ (3/2)·d·log(l)`**

★★★これが**残っているすべて**である。

## ★機構 —— `§9-1023` の共体積表示から

`12·d·ht^Falt(E) = −12·Σᶠ_p neronExp_p·log N(p) − 12·d·log(2π)
                    + 12·Σ_σ log‖u_σ‖ − 6·Σ_σ log(covol P_σ)`（`§9-1023`、**無条件**）

★差を取ると `log(2π)` の項が消え、3 項が残る:

| 項 | 状態 |
|---|---|
| `covol` の比 | ★★**済**——解析的な同種写像で指数 `l` なら `1/l` 倍（`§9-1017`・`§9-1019`） |
| `neronExp` の差 | ☆残る（`hloc` の第 2 項） |
| `u` の比 | ☆残る（`hloc` の第 1 項） |

★★残る 2 項がちょうど `[DelSB616] Théorème 2.4` の**段 1・段 2**
（`ω_{E′} ⊆ ω_E` で余核が `#H` で消える）である。

## ★★入力 `hcov` について

`hcov : ∀ σ, covol (P σ) = l · covol (P′ σ)` は「解析的な同種写像の指数が `l`」の
言い換えであり、`Found/GenEll/IsogenyPeriodPair.lean` の
`covol_eq_index_mul_pair`（`§9-1019`、**無条件**）が与える形である。

## ☆これは仮定を受けているのではない

☆`hloc` は**受けた仮定ではなく、残っている定理の statement そのもの**である。
★本ファイルの意味は「残りがこれ **1 本** であることを型で示す」ことにある
——`.src` は条つきにして指標には数えない（`check.mjs` B6）。
-/

namespace ABC3.Found.GaloisRep

open NumberField IsDedekindDomain WeierstrassCurve ABC3.Found.GenEll
open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-! ## ★★★★★★★★★★★★★★★★★★★★★★帰着 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★**同種写像評価は 1 つの局所不等式に帰着する**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`hcov`（解析的な指数が `l`）は `§9-1019` の `covol_eq_index_mul_pair` が与える。
☆`hloc` が**残っているすべて**である——`u` の比と `neronExp` の差。

★★★これが `Found/GaloisRep/Lemma35Concrete.lean` の仮説 `hfalt` を埋める形であり、
`hfalt` が埋まれば `Lemma 3.5` → `Lemma 3.7` → `Theorem 3.8` →
`Corollary 4.3` / `4.4` が続く。 -/
theorem htFalt_isogeny_le (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (l : ℕ) (hl : 1 ≤ l)
    (P P' : (L →+* ℂ) → PeriodPair) (C C' : (L →+* ℂ) → VariableChange ℂ)
    (hPC : ∀ σ, C σ • (E.map σ) = latticeCurve (P σ))
    (hPC' : ∀ σ, C' σ • (E'.map σ) = latticeCurve (P' σ))
    (hcov : ∀ σ, covol (P σ) = (l : ℝ) * covol (P' σ))
    (hloc : ((∑ σ : (L →+* ℂ), Real.log ‖((C' σ).u : ℂ)‖)
              - (∑ σ : (L →+* ℂ), Real.log ‖((C σ).u : ℂ)‖))
            - ((∑ᶠ p : HeightOneSpectrum (𝓞 L),
                  (neronExp p E' : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
               - (∑ᶠ p : HeightOneSpectrum (𝓞 L),
                  (neronExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal)))
          ≤ (3/2) * (Module.finrank ℚ L : ℝ) * Real.log l) :
    htFaltOf L E' ≤ htFaltOf L E + 2 * Real.log l := by
  have hd : (0:ℝ) < (Module.finrank ℚ L : ℝ) := by exact_mod_cast Module.finrank_pos
  have hl0 : (0:ℝ) < (l:ℝ) := by exact_mod_cast hl
  have hrep := twelve_finrank_htFaltOf_eq E P C hPC
  have hrep' := twelve_finrank_htFaltOf_eq E' P' C' hPC'
  have hV : (∑ σ : (L →+* ℂ), Real.log (covol (P σ)))
      = (Module.finrank ℚ L : ℝ) * Real.log l
        + ∑ σ : (L →+* ℂ), Real.log (covol (P' σ)) := by
    have hterm : ∀ σ : (L →+* ℂ), Real.log (covol (P σ))
        = Real.log l + Real.log (covol (P' σ)) := by
      intro σ
      rw [hcov σ, Real.log_mul (ne_of_gt hl0) (ne_of_gt (covol_pos (P' σ)))]
    rw [Finset.sum_congr rfl (fun σ _ => hterm σ), Finset.sum_add_distrib,
      Finset.sum_const, Finset.card_univ, NumberField.Embeddings.card L ℂ, nsmul_eq_mul]
  rw [hV] at hrep
  have hkey : 12 * (Module.finrank ℚ L : ℝ) * htFaltOf L E'
      ≤ 12 * (Module.finrank ℚ L : ℝ) * (htFaltOf L E + 2 * Real.log l) := by
    have hexp : 12 * (Module.finrank ℚ L : ℝ) * (htFaltOf L E + 2 * Real.log l)
        = 12 * (Module.finrank ℚ L : ℝ) * htFaltOf L E
          + 24 * ((Module.finrank ℚ L : ℝ) * Real.log l) := by ring
    rw [hrep', hexp, hrep]
    linarith [hloc]
  have h12 : (0:ℝ) < 12 * (Module.finrank ℚ L : ℝ) := by linarith
  nlinarith [hkey, h12]

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★正準形——選択に依らない形 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**残り 1 本の正準形**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

    **`(l−1)·d·deg∞(E) − (archSum(E′) − archSum(E)) ≤ 24·d·log(l)`**

★上の `htFalt_isogeny_le` は周期対の選択 `(P_σ, C_σ)` に依る形だったが、
**こちらは本プロジェクトが定義した量 `degInfOf`・`archSum` だけで書ける**
——選択に依らない。

★★★`§9-1025`（第 581）の測定どおり、`htFalt_isogeny_le` の `hloc` は
`§9-1021`〜`§9-1023` の道具を通すとこの形に**戻る**（一周する）。
したがって**これが残り 1 本の正準形**である。

☆`hdeg` は `Lemma 3.2, (ii)`（`§9-1011`）＋その大域化（`§9-1012`）が与える。
☆`hArch` が残っているすべてであり、`[DelSB616] Théorème 2.4` の段 1・段 2 に当たる。 -/
theorem htFalt_isogeny_le_canonical (E E' : WeierstrassCurve L) (l : ℕ)
    (hdeg : degInfOf L E' = (l : ℝ) * degInfOf L E)
    (hArch : ((l : ℝ) - 1) * (Module.finrank ℚ L : ℝ) * degInfOf L E
        - (archSum L E' - archSum L E)
      ≤ 24 * (Module.finrank ℚ L : ℝ) * Real.log l) :
    htFaltOf L E' ≤ htFaltOf L E + 2 * Real.log l := by
  have hd : (0:ℝ) < (Module.finrank ℚ L : ℝ) := by exact_mod_cast Module.finrank_pos
  have h1 : 12 * htFaltOf L E' = degInfOf L E' - archSum L E' / (Module.finrank ℚ L : ℝ) := by
    rw [htFaltOf]; field_simp
  have h2 : 12 * htFaltOf L E = degInfOf L E - archSum L E / (Module.finrank ℚ L : ℝ) := by
    rw [htFaltOf]; field_simp
  have hdiff : 12 * htFaltOf L E' - 12 * htFaltOf L E
      = ((l : ℝ) - 1) * degInfOf L E
        - (archSum L E' - archSum L E) / (Module.finrank ℚ L : ℝ) := by
    rw [h1, h2, hdeg]; ring
  have hmul : (((l : ℝ) - 1) * degInfOf L E
      - (archSum L E' - archSum L E) / (Module.finrank ℚ L : ℝ))
      * (Module.finrank ℚ L : ℝ)
      = ((l : ℝ) - 1) * (Module.finrank ℚ L : ℝ) * degInfOf L E
        - (archSum L E' - archSum L E) := by
    field_simp
  have hkey : ((l : ℝ) - 1) * degInfOf L E
      - (archSum L E' - archSum L E) / (Module.finrank ℚ L : ℝ) ≤ 24 * Real.log l := by
    refine le_of_mul_le_mul_right ?_ hd
    rw [hmul]
    nlinarith [hArch]
  linarith [hdiff, hkey]

/-! ## ★出典の紐付け(`.src`)——★★**条つき（`hloc`／`hArch` が残っている）** -/

def htFalt_isogeny_le_canonical.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(残り 1 本の正準形——deg∞ と archSum だけで書ける)",
    sectionId := "genell-prop-3-4" }

def htFalt_isogeny_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(同種写像評価は 1 つの局所不等式に帰着する)",
    sectionId := "genell-prop-3-4" }

def htFalt_isogeny_le.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "twelve_finrank_htFaltOf_eq(ht^Falt の共体積表示、§9-1023)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.twelve_finrank_htFaltOf_eq") 4,
    .citation "[ABC3]" "covol_eq_index_mul_pair(hcov を与える、§9-1019)"
      (.inProject "ABC3" "ABC3.Found.GenEll.covol_eq_index_mul_pair") 3,
    .otherPaper "[DelSB616]"
      ("Théorème 2.4 の段 1・段 2 —— ω_{E′} ⊆ ω_E で余核が #H で消えること。" ++
       "☆★★★★★★★★★★**hloc がその内容であり、§3・§4 の残り 5 項目すべてが" ++
       "これ 1 本を待っている**。" ++
       "★Skeleton/GenEll/IsogenyHeight.lean に節点があり、" ++
       "0_Source に原典がある(2026-08-29 の発見)") 14,
    .implicitStep
      ("★★★★★★★★★★到達点(2026-08-29、第 580): " ++
       "§9-1023 の共体積表示から差を取ると log(2π) の項が消え、3 項が残る" ++
       "——covol の比・neronExp の差・u の比。" ++
       "★**covol の比は取れている**(§9-1017・§9-1019、指数 l で 1/l 倍)。" ++
       "★★したがって残るのは hloc ただ 1 本であり、" ++
       "それがちょうど [DelSB616] の段 1・段 2 である") 9,
    .implicitStep
      ("☆hloc は**受けた仮定ではなく、残っている定理の statement そのもの**である。" ++
       "★本ファイルの意味は「残りがこれ 1 本であることを型で示す」ことにあり、" ++
       ".src は条つきにして指標には数えない(check.mjs B6)") 8 ]

end ABC3.Found.GaloisRep
