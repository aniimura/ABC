import ABC3.Meta.Claim
import ABC3.Skeleton.GenEll.Section2
import ABC3.Skeleton.GenEll.Section4
import ABC3.Skeleton.IUTchIII.Cor312
import ABC3.Skeleton.PGC.Section4
import Mathlib.Data.Nat.Factorization.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real

/-!
# 最終目標 —— abc 予想と、各論文の頂点を繋ぐ 1 本の鎖

★このファイルの目的は **証明ではなく接続**である。

これまで `Skeleton/` は 16 本の論文の主張が**それぞれ島**として型になっている状態だった
(実測 2026-09-09: `theorem_4_2` を消費する宣言は 0、abc にあたる宣言は木に 0 件)。
★**繋ぐ先が無ければ「型が最後まで繋がっているか」は問えない。**

そこで本ファイルは 3 つを置く。

1. **最終目標**を、パラメータを持たない閉じた命題として書く(`ABCConjecture`)。
2. **各論文の頂点**を、その論文の `theorem` の型**そのもの**として `Prop` に名付ける
   (`Top_*`)。★`:= @<木の定理>` が型検査を通ることが、
   「この命題が木の主張と字面まで一致する」ことの**機械的証拠**である
   ——写し間違えればビルドが落ちる。
3. 頂点から目標への**還元 1 本**(`abc_of_tops`)。★これは `sorry` である。

★★**`sorry` の位置がそのまま「不足」の位置である。**
`#print axioms abcConjecture_holds` と `node tools/goal-chain.mjs` で、
目標が現在どの未証明宣言に載っているかが機械で列挙できる。

## ★このファイルが**主張していないこと**

- `abc_of_tops` の数学的内容を本体は**読んでいない**。原典の該当箇所
  ([GenEll] Theorem 2.1 の (ii)⟹(i) 側 + [IUTchIV] Theorem 1.10)を
  写したものではなく、**空欄**である。
- `Top_*` が**空虚でない**ことは、このファイルからは分からない。
  `∀ D : EllModuliData, …` は `EllModuliData` が 1 つも構成されていなければ
  自明に真になりうる。★その検査は `tools/goal-chain.mjs --vacuity` が別に行う。
-/

namespace ABC3.Skeleton.Goal

open ABC3.Meta

/-! ## ① 最終目標 -/

/-- 根基 `rad n = ∏_{p ∣ n} p`。 -/
def rad (n : ℕ) : ℕ := ∏ q ∈ n.primeFactors, q

/-- **abc 予想**(Masser–Oesterlé)。

> 任意の `ε > 0` に対し定数 `C` が存在して、互いに素な正整数 `a + b = c` について
> `c ≤ C · rad(abc)^{1+ε}`。

★パラメータを持たない閉じた命題である——`Interface` の欄も `Data` 構造体も通らない。
★したがって**空虚になりようがない**。木で唯一この性質を持つ命題である。 -/
def ABCConjecture : Prop :=
  ∀ ε : ℝ, 0 < ε → ∃ C : ℝ, 0 < C ∧
    ∀ a b c : ℕ, 0 < a → 0 < b → a + b = c → Nat.Coprime a b →
      (c : ℝ) ≤ C * (rad (a * b * c) : ℝ) ^ (1 + ε)

/-! ## ② 各論文の頂点 —— 型は木の定理そのもの

★各 `Top_*` の直後の `theorem top_* : Top_* := @<木の定理>` が**接続の検査**である。 -/

section Tops

open ABC3.Interface.GenEll ABC3.Interface.IUTchIII ABC3.Interface.PGC
open ABC3.Found.GenEll

/-- **[GenEll] Theorem 2.1** の型。 -/
def Top_GenEll_Thm21 : Prop :=
  ∀ {P : Type} (T : Thm21Data P) (htOmega logDiff logCond : P → ℝ) (eps : ℝ),
    BDle (fun x : ↑T.degLe => (1 + eps) * (logDiff x.1 + logCond x.1))
         (fun x : ↑T.degLe => htOmega x.1)
    ↔
    (∀ KV : Set P, T.CB KV →
      BDle (fun x : ↑(KV ∩ T.degLe) => (1 + eps) * (logDiff x.1 + logCond x.1))
           (fun x : ↑(KV ∩ T.degLe) => htOmega x.1))

theorem top_genEll_thm21 : Top_GenEll_Thm21 := @ABC3.Skeleton.GenEll.theorem_2_1

/-- **[GenEll] Corollary 4.4** の型(★GenEll トラックの北極星)。 -/
def Top_GenEll_Cor44 : Prop :=
  ∀ (D : EllModuliData) (KV : Set D.EllClass), D.CompactlyBounded KV →
    ∃ C : ℝ, 0 < C ∧ ∃ Exc : Set D.EllClass, D.GaloisFinite Exc ∧
      ∀ (E : D.Curve) (S : Finset ℕ), (∀ p ∈ S, p.Prime) →
        D.MinimalField E → D.cls E ∈ KV → D.cls E ∉ Exc →
        ∃ lo lb : ℕ, Nat.Prime lo ∧ Nat.Prime lb ∧ lo ∉ S ∧ lb ∉ S
          ∧ (D.PrimeToMultPrimes E lo ∧ D.PrimeToLocalHeights E lo
              ∧ D.PrimeToMultPrimes E lb ∧ D.PrimeToLocalHeights E lb
              ∧ D.PrimeToRamification E lb)
          ∧ (D.ImageContainsSL2 E lo ∧ D.ImageSurjective E lb)
          ∧ (lo : ℝ) ≤ 23040 * 100 * (D.degOfDefinition E : ℝ) * D.faltingsHeight (D.cls E)
              + 2 * (∑ p ∈ S, Real.log p) + C * (D.degOfDefinition E : ℝ)
          ∧ (lb : ℝ) ≤ 23040 * 100 * (D.degOfDefinition E : ℝ) * D.faltingsHeight (D.cls E)
              + 6 * (D.degOfDefinition E : ℝ) * D.logDiffMell (D.cls E)
              + 2 * (∑ p ∈ S, Real.log p) + C * (D.degOfDefinition E : ℝ)

theorem top_genEll_cor44 : Top_GenEll_Cor44 := @ABC3.Skeleton.GenEll.cor_4_4

/-- **[IUTchIII] Corollary 3.12** の型。 -/
def Top_IUT_Cor312 : Prop :=
  ∀ (D : PilotObjectData),
    ABC3.Skeleton.IUTchIII.PossibleImagesContained D →
    ABC3.Skeleton.IUTchIII.LogShellPacketCompact D →
    ABC3.Skeleton.IUTchIII.HullCompactOfRelCompact D →
    ABC3.Skeleton.IUTchIII.OutputLogVolumesEq D →
    ABC3.Skeleton.IUTchIII.QLogVolMem D →
    ABC3.Skeleton.IUTchIII.thetaLogVol D ≠ ⊤
      ∧ ABC3.Skeleton.IUTchIII.qLogVol D ≤ ABC3.Skeleton.IUTchIII.thetaLogVol D

theorem top_iut_cor312 : Top_IUT_Cor312 := @ABC3.Skeleton.IUTchIII.cor_3_12

/-- **[pGC] Theorem 4.2** の型。★木では `sorry`(全射性が残っている)。 -/
def Top_PGC_Thm42 : Prop :=
  ∀ {p : ℕ} [Fact p.Prime] (K K' : ABC3.Skeleton.PGC.PAdicLocalField p),
    Function.Bijective
      (ABC3.Found.PGC.naturalOuterIso (ABC3.Found.PGC.ramificationFiltration p)
        (ABC3.Found.PGC.isNaturalFiltration_ramificationFiltration p) (K := K) (K' := K'))

theorem top_pgc_thm42 : Top_PGC_Thm42 := @ABC3.Skeleton.PGC.theorem_4_2

end Tops

/-! ## ③ 還元 —— ★ここが空欄である -/

/-- **頂点 4 本から abc へ**。

★★**本体はこの還元の証明を読んでいない。`sorry` である。**

原典側で対応する道筋(★未検証、`.needs` に読むべき箇所だけを書く):
[IUTchIII] Cor 3.12 → [IUTchIV] Theorem 1.10(高さ不等式)→
[GenEll] Corollary 4.4(compactly bounded の場合)→
[GenEll] Theorem 2.1 の `(ii) ⟹ (i)`(一般の場合へ)→ `X = ℙ¹ ∖ {0,1,∞}`、`d = 1` で abc。

★[pGC] Theorem 4.2 は上の鎖の**どこで使われるか**を本体はまだ特定していない
(遠アーベル幾何の入力として [AbsTopIII] 経由で入るはずだが、**測っていない**)。
それでも仮説として並べてあるのは、★**「使われない頂点」を後で機械が検出できるようにする**ためである。 -/
theorem abc_of_tops
    (_h21 : Top_GenEll_Thm21) (_h44 : Top_GenEll_Cor44)
    (_h312 : Top_IUT_Cor312) (_h42 : Top_PGC_Thm42) :
    ABCConjecture := sorry

/-- ★**最終目標の宣言**。

`#print axioms abcConjecture_holds` が `sorryAx` を出す間、abc は未証明である。
★どの `sorry` に載っているかは `node tools/goal-chain.mjs` が列挙する。 -/
theorem abcConjecture_holds : ABCConjecture :=
  abc_of_tops top_genEll_thm21 top_genEll_cor44 top_iut_cor312 top_pgc_thm42

/-! ## ★非空虚性の対照(G9)

`ABCConjecture` は**条件を 1 つも持たない閉じた命題**なので空虚になりようがない。
それでも内側の全称が空でないことを 1 例で示す —— `1 + 8 = 9`、`gcd(1,8) = 1`。 -/

theorem abcConjecture_holds.nonvacuous (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, 0 < C ∧ (9 : ℝ) ≤ C * (rad (1 * 8 * 9) : ℝ) ^ (1 + eps) := by
  obtain ⟨C, hC, h⟩ := abcConjecture_holds eps heps
  exact ⟨C, hC, h 1 8 9 (by norm_num) (by norm_num) (by norm_num) (by decide)⟩

/-- `abc_of_tops` の 4 つの仮説がどれも空でないことの対照 —— 実際に木の頂点で埋まる。 -/
theorem abc_of_tops.nonvacuous : ABCConjecture :=
  abc_of_tops top_genEll_thm21 top_genEll_cor44 top_iut_cor312 top_pgc_thm42

/-! ## 出典 -/

def rad.src : Source :=
  { paper := "GenEll", pdfPage := 11, item := "Theorem 2.1", sectionId := "genell-thm-2-1" }

def ABCConjecture.src : Source :=
  { paper := "GenEll", pdfPage := 11, item := "Theorem 2.1", sectionId := "genell-thm-2-1" }

def Top_GenEll_Thm21.src : Source :=
  { paper := "GenEll", pdfPage := 11, item := "Theorem 2.1 (頂点の型としての引用)", sectionId := "genell-thm-2-1" }

def top_genEll_thm21.src : Source :=
  { paper := "GenEll", pdfPage := 11, item := "Theorem 2.1 (頂点の型としての引用)", sectionId := "genell-thm-2-1" }

def Top_GenEll_Cor44.src : Source :=
  { paper := "GenEll", pdfPage := 23, item := "Corollary 4.4 (頂点の型としての引用)", sectionId := "genell-cor-4-4" }

def top_genEll_cor44.src : Source :=
  { paper := "GenEll", pdfPage := 23, item := "Corollary 4.4 (頂点の型としての引用)", sectionId := "genell-cor-4-4" }

def Top_IUT_Cor312.src : Source :=
  { paper := "IUTchIII", pdfPage := 174, item := "Corollary 3.12 (頂点の型としての引用)",
    sectionId := "cor-3-12-conclusion" }

def top_iut_cor312.src : Source :=
  { paper := "IUTchIII", pdfPage := 174, item := "Corollary 3.12 (頂点の型としての引用)",
    sectionId := "cor-3-12-conclusion" }

def Top_PGC_Thm42.src : Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2 (頂点の型としての引用)", sectionId := "theorem-4-2" }

def top_pgc_thm42.src : Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2 (頂点の型としての引用)", sectionId := "theorem-4-2" }

def abcConjecture_holds.src : Source :=
  { paper := "GenEll", pdfPage := 11, item := "Theorem 2.1", sectionId := "genell-thm-2-1" }

def Top_GenEll_Thm21.needs : List ProofObligation := []
def top_genEll_thm21.needs : List ProofObligation := []
def Top_GenEll_Cor44.needs : List ProofObligation := []
def top_genEll_cor44.needs : List ProofObligation := []
def Top_IUT_Cor312.needs : List ProofObligation := []
def top_iut_cor312.needs : List ProofObligation := []
def Top_PGC_Thm42.needs : List ProofObligation := []
def top_pgc_thm42.needs : List ProofObligation := []

def abcConjecture_holds.needs : List ProofObligation :=
  [ .implicitStep
      ("★abc_of_tops と同じ 1 点のみ。頂点 4 本は型検査で接続済みである" ++
       "(`:= @<木の定理>` が通ることがその証拠)。") 11 ]

def abc_of_tops.src : Source :=
  { paper := "GenEll", pdfPage := 11, item := "Theorem 2.1", sectionId := "genell-thm-2-1" }

def abc_of_tops.needs : List ProofObligation :=
  [ .implicitStep
      ("★本体が読んでいない段。[IUTchIV] Theorem 1.10 から [GenEll] Corollary 4.4 への" ++
       "翻訳、および Theorem 2.1 の (ii) ⟹ (i) を X = ℙ¹∖{0,1,∞}, d = 1 に当てる段。" ++
       "★工数は測っていない。") 11
  , .implicitStep
      ("★[pGC] Theorem 4.2 が鎖のどこに入るかの特定。仮説には並べたが、" ++
       "本体は接続点を測っていない。") 7 ]

end ABC3.Skeleton.Goal
