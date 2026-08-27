/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.BDClass

/-!
# [GenEll] Proposition 3.4 の有限性の主張 —— **不等式と Northcott から**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> on Mell(Q). In particular, if C ∈R, then the set of points [E] ∈Mell(Q)≤d such

## ★★★原文が「In particular」で畳んだ段

`Proposition 3.4` は 2 つを言う:

1. `deg ≾ ht_∞ ≾ 12(1+ε)·ht_Falt ≾ (1+ε)·ht_∞`（★`[Silv2]` に依る——未）
2. **「In particular」** —— `ht_Falt ≤ C` なる点は有限個

★★**2 は 1 と Northcott から初等的に出る**。本ファイルはその段を取る。

## ★★★★機構は 2 行

`ht_∞ ≾ a·ht_Falt`（`a ≥ 0`、定数差込み）なら

    ht_Falt p ≤ C  ⟹  a·ht_Falt p ≤ a·C  ⟹  ht_∞ p ≤ a·C + C′

なので `{ht_Falt ≤ C} ⊆ {ht_∞ ≤ a·C + C′}`。★あとは `ht_∞` の Northcott 性である。

## ★★★★★2 つの入力がどこにあるか

| 入力 | 状態 |
|---|---|
| `ht_∞ ≾ a·ht_Falt` | ★`Proposition 3.4` の 1（`[Silv2]` Prop 2.1 に依る——未） |
| `ht_∞` の Northcott 性 | ★★`Proposition 1.4, (iv)`。`M_ell` の粗モジュライは `j`-線なので `ℙ¹` の場合で足りる（`Check/GenEll/NorthcottProjModelNonvacuous.lean` で検査済み） |

★★★**したがって残るのは `[Silv2]` の側だけ**である
（`Found/GenEll/Covolume.lean` ほかで並行セッションが構築中）。

## ★逸脱の記録（CLAUDE.md の「逸脱」）

★本ファイルは**点の型 `P` と実数値関数**で述べる——`M_ell(ℚ̄)^{≤d}` という
具体的な型は取らない。★★原文の主張をこの形に写したのは、
`Interface` の `EllModuli` に依らずに段を取るためである。
-/

namespace ABC3.Found.GenEll

/-- ★★★★★★**[GenEll] Proposition 3.4 の有限性の主張**（原文の `≾` の形）。

原文 (GenEll p.17):
> on Mell(Q). In particular, if C ∈R, then the set of points [E] ∈Mell(Q)≤d such

★`ht_∞ ≾ a·ht_Falt` と `ht_∞` の Northcott 性から、`ht_Falt` の有限性が出る。

★★機構は包含 1 本: `{ht_Falt ≤ C} ⊆ {ht_∞ ≤ a·C + C′}`。 -/
theorem finite_of_le_of_northcott {P : Type*} (htFalt htInf : P → ℝ)
    (a : ℝ) (ha : 0 ≤ a) (C' : ℝ)
    (hle : ∀ p, htInf p ≤ a * htFalt p + C')
    (hN : ∀ B : ℝ, {p | htInf p ≤ B}.Finite)
    (C : ℝ) : {p | htFalt p ≤ C}.Finite := by
  refine Set.Finite.subset (hN (a * C + C')) (fun p hp => ?_)
  have h1 : htFalt p ≤ C := hp
  have h2 : a * htFalt p ≤ a * C := mul_le_mul_of_nonneg_left h1 ha
  exact le_trans (hle p) (by linarith)

/-- ★★★★★**下からの形** —— `c·ht_∞ − C′ ≤ ht_Falt` からも同じ結論。

★原文の鎖 `ht_∞ ≾ 12(1+ε)·ht_Falt` はこちらの形でも読める（`c = 1/(12(1+ε))`）。 -/
theorem finite_of_bdge_of_northcott {P : Type*} (htFalt htInf : P → ℝ)
    (c : ℝ) (hc : 0 < c) (C' : ℝ)
    (hle : ∀ p, c * htInf p - C' ≤ htFalt p)
    (hN : ∀ B : ℝ, {p | htInf p ≤ B}.Finite)
    (C : ℝ) : {p | htFalt p ≤ C}.Finite := by
  refine Set.Finite.subset (hN ((C + C') / c)) (fun p hp => ?_)
  have h1 : c * htInf p - C' ≤ C := le_trans (hle p) hp
  show htInf p ≤ (C + C') / c
  rw [le_div_iff₀ hc]
  linarith

/-- ★★**`BDle` から不等式を取り出す** —— 原文の `≾` を上の補題の仮定の形に直す。

★`BDle α β` は `∃ C, ∀ x, β x − α x ≤ C` である（`BDClass.lean`）——**向きに注意**。 -/
theorem exists_le_of_bdle {P : Type*} (f g : P → ℝ) (h : BDle f g) :
    ∃ C' : ℝ, ∀ p, g p ≤ f p + C' := by
  obtain ⟨C', hC'⟩ := h
  exact ⟨C', fun p => by linarith [hC' p]⟩

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` は置かない。** `Proposition 3.4` の主張 1
（`deg ≾ ht_∞ ≾ 12(1+ε)·ht_Falt ≾ (1+ε)·ht_∞`）は `[Silv2]` Prop 2.1 に依り、未である。 -/

def finite_of_le_of_northcott.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(有限性の主張——不等式と Northcott から。不等式自体は含まない)",
    sectionId := "genell-prop-3-4" }

def finite_of_le_of_northcott.needs : List ABC3.Meta.ProofObligation :=
  [ .otherPaper "[GenEll]"
      ("Proposition 1.4, (iv)(Northcott 性)——M_ell の粗モジュライは j-線なので " ++
       "ℙ¹ の場合で足りる。Check/GenEll/NorthcottProjModelNonvacuous.lean で検査済み") 6,
    .citation "[Silv2]" "Proposition 2.1(deg ≾ ht_∞ ≾ 12(1+ε)·ht_Falt の側)"
      (.absent "0_Source に [Silv2] は無く、mathlib にも Faltings 高さと無限素点での高さの比較は無い") 17,
    .implicitStep
      ("★原文が「In particular」で畳んだ段である。機構は包含 1 本: " ++
       "{ht_Falt ≤ C} ⊆ {ht_∞ ≤ a·C + C′}") 17,
    .implicitStep
      ("★逸脱: 点の型 P と実数値関数で述べる——M_ell(ℚ̄)^{≤d} という具体的な型は取らない。" ++
       "Interface の EllModuli に依らずに段を取るためである") 17 ]

def finite_of_bdge_of_northcott.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(有限性の主張——下からの形)",
    sectionId := "genell-prop-3-4" }

def exists_le_of_bdle.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(BDle から不等式を取り出す)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
