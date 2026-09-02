/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.Lemma35Assemble
import ABC3.Found.GenEll.LCyclicPoint
import ABC3.Found.GenEll.Lemma37StableLineCop
import ABC3.Meta.Claim

/-!
# 第 1225 ブロック —— **`hdag` を `DegCurve` の語彙で出す**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17–p.18。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か

第 1223（`lemma_3_5_height_ineq_stableLine`）を `DegCurve` の語彙に読み替え、
第 1210（`lemma_3_7_stableLine_cop`）が要る `hdag` の形にする。

| 段 | 材料 | 第 |
|---|---|---|
| `HasLCyclicJ` → 安定な位数 `l` の点 | `exists_stablePoint_of_hasLCyclicJ` | 1205 |
| `PrimeToLocalHeights` → `l ∤ v_p(j)` | `not_dvd_jExp_of_primeToLocalHeights` | 1151 |
| 不等式 | `lemma_3_5_height_ineq_stableLine` | 1223 |
| `1/(12(1+ε))` → `1/14` | `ε ≔ 1/6` | — |

★★逸脱: Vélu の商が楕円かつ半安定であることは `VeluQuotOK` として受ける。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Interface.GaloisRep WeierstrassCurve
open IsDedekindDomain NumberField Finset ABC3.Meta
open scoped Classical

/-- ☆Vélu の商が楕円かつ半安定であること（原文が「同種なので自動」と括弧で述べる段）。

★★★★★★★★☆**2026-09-02（第 1322-1323）——悪い素点側の見通し**

☆半安定性は 2 つに分かれる:

* **悪い素点**（`jExp < 0`）——★**核は取れた**。
  `isUnit_c4_velu_tate`（第 1323）が「Tate 曲線の Vélu の商は `c₄` が単元」を、
  `semistableAt_of_c4_valAdd_zero`（第 1322）が「`c₄` が単元の整モデルは半安定」を与える。
  ☆残るのは商のモデルを Tate モデルに移す変数変換
  （`veluQuotientFull_vcPoint_eq` 第 969・`vAdd_tateModel_u_eq_zero` 第 1056、在庫）と
  整性（`veluIntegralClosed`、在庫）の配管である。
* **良い素点**（`jExp ≥ 0`）——☆**同種で良還元が保たれる**（Néron–Ogg–Shafarevich）が要る。

★これと Vélu の商の楕円性（Vélu の定理）が `VeluQuotOK` に残る 2 つの既知数学である。 -/
def VeluQuotOK (E : SSCurve) (l : ℕ) : Prop :=
  ∀ (M : IntermediateField E.fld E.alg) [FiniteDimensional E.fld M]
     (Q' : (E.W.baseChange M).toAffine.Point),
     letI : DecidableEq (M : Type) := fun a b => Classical.propDecidable (a = b)
     letI : NumberField (M : Type) := NumberField.of_module_finite E.fld M
     letI : IsScalarTower (𝓞 E.fld) E.fld M := isScalarTower_ringOfIntegers_base E.fld M
     letI : IsScalarTower (𝓞 E.fld) (𝓞 (M : Type)) M := isScalarTower_ringOfIntegers_top E.fld M
     (veluQuotientFull (E.W.baseChange M)
        (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q')))).IsElliptic ∧
     ∀ P : HeightOneSpectrum (𝓞 (M : Type)),
       SemistableAt P (veluQuotientFull (E.W.baseChange M)
         (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q'))))

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**`hdag` を `DegCurve` の語彙で出す**（第 1225）。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

☆`HasLCyclicJ` と `PrimeToLocalHeights` と `VeluQuotOK` から
`Lemma 3.5` の結論 `(†)` が出る。

★★★これが第 1210（`lemma_3_7_stableLine_cop`）の `hdag` そのものである。 -/
theorem hdag_of_stableLine :
    ∃ C₅ : ℝ, 0 ≤ C₅ ∧
      ∀ (E : DegCurve) (l : ℕ), Nat.Prime l → HasLCyclicJ E.toSSCurve l →
        E.toSSCurve.PrimeToLocalHeights l → VeluQuotOK E.toSSCurve l →
        ((l : ℝ) / 14) * degInfOf E.toSSCurve.fld E.toSSCurve.W
          ≤ htFaltOf E.toSSCurve.fld E.toSSCurve.W + 2 * Real.log l + C₅ := by
  obtain ⟨C, hC⟩ := lemma_3_5_height_ineq_stableLine (1 / 6) (by norm_num)
  refine ⟨max C 0, le_max_right _ _, ?_⟩
  intro E l hl hcyc hcop hquot
  haveI hEell : ((E.toSSCurve.W).baseChange E.toSSCurve.alg).IsElliptic := by
    show ((E.toSSCurve.W).map (algebraMap E.toSSCurve.fld E.toSSCurve.alg)).IsElliptic
    infer_instance
  obtain ⟨Q, hQ, hst⟩ := exists_stablePoint_of_hasLCyclicJ E.toSSCurve l hl hcyc
  have hmain := hC E.toSSCurve.fld E.toSSCurve.W l hl E.toSSCurve.ss
    (E.toSSCurve.not_dvd_jExp_of_primeToLocalHeights hl hcop) Q hQ hst hquot
  have hcoef : (1 : ℝ) / (12 * (1 + 1 / 6)) = 1 / 14 := by norm_num
  rw [hcoef] at hmain
  have hle : C ≤ max C 0 := le_max_left _ _
  nlinarith [hmain, hle]

/-! ## ★★★★★★★★★★★★★★★★★★★★`Lemma 3.7` を安定直線の側で閉じる -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] Lemma 3.7 —— 安定直線の側（`hdag` を埋めた形）**（第 1226）。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

☆第 1210 の `hdag` を第 1225 で埋めた。
★受けている外部の仮定は **`VeluQuotOK` ただ 1 つ**である
——Vélu の商が楕円かつ半安定であること（原文が「同種なので自動」と括弧で述べる段）。

★★★これで `Theorem 3.8` が要る `¬ HasLCyclicJ` が `[E] ∉ Exc` から出る。 -/
theorem lemma_3_7_stableLine_full (KV : Set ℂ) (hKV : CompactlyBoundedJ KV)
    (eps : ℝ) (heps : 0 < eps)
    (hquot : ∀ (E : SSCurve) (l : ℕ), VeluQuotOK E l) :
    ∃ C : ℝ, 0 < C ∧ ∃ Exc : Set ℂ, GaloisFiniteJ Exc ∧
      ∀ (E : DegCurve) (l : ℕ), Nat.Prime l →
        ∀ condA condB : Prop,
          (condA ↔ (100 * (E.deg : ℝ)
                      * (faltingsHeightJ E.j + C * (E.deg : ℝ) ^ eps) ≤ (l : ℝ)
                    ∧ E.toSSCurve.HasMultRed)) →
          (condB ↔ (E.j ∈ KV ∧ E.toSSCurve.PrimeToLocalHeights l)) →
          (condA → E.toSSCurve.PrimeToLocalHeights l)
        ∧ (condB → E.j ∉ Exc → E.toSSCurve.HasMultRed)
        ∧ ((condA ∨ condB) → HasLCyclicJ E.toSSCurve l → E.j ∈ Exc) := by
  obtain ⟨C₅, hC₅, hdag⟩ := hdag_of_stableLine
  exact lemma_3_7_stableLine_cop KV hKV eps heps C₅ hC₅
    (fun E l hl hcyc hcop => hdag E l hl hcyc hcop (hquot _ _))

/-! ## ★★★★★★★★★★★★`lcyclicExc` の形 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★
**`l`-巡回な `j` は Galois-有限な例外集合に入る**（第 1239）。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

☆第 1226 の第 3 の主張を取り出した形。
★★★これが `EllModuliWitness` の `lcyclicExc` 欄
（`galoisFiniteJ_lcyclicExcJ`）が要る形である。 -/
theorem exists_galoisFinite_lcyclic (KV : Set ℂ) (hKV : CompactlyBoundedJ KV)
    (eps : ℝ) (heps : 0 < eps)
    (hquot : ∀ (E : SSCurve) (l : ℕ), VeluQuotOK E l) :
    ∃ C : ℝ, 0 < C ∧ ∃ Exc : Set ℂ, GaloisFiniteJ Exc ∧
      ∀ (E : DegCurve) (l : ℕ), Nat.Prime l →
        ((100 * (E.deg : ℝ)
              * (faltingsHeightJ E.j + C * (E.deg : ℝ) ^ eps) ≤ (l : ℝ)
            ∧ E.toSSCurve.HasMultRed)
          ∨ (E.j ∈ KV ∧ E.toSSCurve.PrimeToLocalHeights l)) →
        HasLCyclicJ E.toSSCurve l → E.j ∈ Exc := by
  obtain ⟨C, hCpos, Exc, hExc, h37⟩ := lemma_3_7_stableLine_full KV hKV eps heps hquot
  refine ⟨C, hCpos, Exc, hExc, ?_⟩
  intro E l hl hcond hcyc
  exact (h37 E l hl _ _ Iff.rfl Iff.rfl).2.2 hcond hcyc

/-! ## ★出典の紐付け(`.src`) -/

def exists_galoisFinite_lcyclic.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(l-巡回な j は Galois-有限な例外集合に入る——lcyclicExc 欄の形)",
    sectionId := "genell-lemma-3-7" }

def lemma_3_7_stableLine_full.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(安定直線の側——hdag を埋めた形。外部の仮定は VeluQuotOK のみ)",
    sectionId := "genell-lemma-3-7" }

def lemma_3_7_stableLine_full.needs : List ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_7_stableLine_cop(第 1210、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.lemma_3_7_stableLine_cop") 1,
    .citation "[ABC3]" "hdag_of_stableLine(第 1225、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.hdag_of_stableLine") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1226）**——受けている外部の仮定は" ++
       "**`VeluQuotOK` ただ 1 つ**である——Vélu の商が楼円かつ半安定であること" ++
       "（原文が「同種なので自動」と括弧で述べる段）。" ++
       "★★★これで `Theorem 3.8` が要る `¬ HasLCyclicJ` が `[E] ∉ Exc` から出る。") 3 ]

def VeluQuotOK.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商が楕円かつ半安定であること——原文が括弧で述べる段)",
    sectionId := "genell-lemma-3-5" }

def hdag_of_stableLine.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(hdag を DegCurve の語彙で出す——商の条件は VeluQuotOK として受ける)",
    sectionId := "genell-lemma-3-7" }

def hdag_of_stableLine.needs : List ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_5_height_ineq_stableLine(第 1223、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.lemma_3_5_height_ineq_stableLine") 1,
    .citation "[ABC3]" "exists_stablePoint_of_hasLCyclicJ(第 1205、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_stablePoint_of_hasLCyclicJ") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1225）**——第 1210（`lemma_3_7_stableLine_cop`）の" ++
       "`hdag` そのものである。☆`1/(12(1+ε))` は `ε ≔ 1/6` で `1/14` になる。" ++
       "★★逸脱: Vélu の商が楕円かつ半安定であることは `VeluQuotOK` として受ける。") 3 ]

end ABC3.Found.GenEll
