/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Sl2Padic
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★`Theorem 3.8` の橋 —— 安定な直線が無ければ非上三角（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.20。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★これは何か

原文 p.20 の最終段は 3 つでできている:

1. `l ∤ v(q)` なら mod `l` 像が `α = (1 1 / 0 1)` を含む（`Lemma 3.2` の直前の局所理論）
2. `l`-巡回部分群スキームを持たなければ mod `l` 像が**非上三角**行列を含む
3. その 2 つから `Lemma 3.1, (iv)` で `SL₂(ℤ_l) ⊆ 像`

★★★`Lemma 3.1, (iv)` は `Sl2Padic.lean` に**実装済み**である。
★★**本ファイルは 2 番の群論の側を取る**——`l`-巡回部分群スキームは
`E[l]` の中の **`Gal`-安定な直線**に対応するので、群論としては

    「安定な直線が無い」⟹「非上三角な行列がある」

だけであり、これは**対偶が自明**である（全部が上三角なら `⟨e₁⟩` が安定）。

## ★★★★★これで `Theorem 3.8` に残るのは

| 段 | 状態 |
|---|---|
| 1 局所理論の行列表示（`α` が像に入る） | ☆残る（Tate 一意化と `Lemma 3.2, (i)` は**閉じている**） |
| ★2 安定な直線が無い ⟹ 非上三角 | ★★**本ファイル**（群論の側）／☆`l`-巡回 ⟷ 安定直線の対応は残る |
| 3 `Lemma 3.1, (iv)` | ★**済み**（`Sl2Padic.lean`） |
| `Lemma 3.7` | ☆`Prop 3.4` 待ち |
| `torsionExt`（3・5 捩れを有理化して半安定に） | ☆残る |

★★`Galois 表現そのもの`は `Found/GaloisRep/GalRep.lean` の `galRep` で**構成済み**である
（2026-08-29 に `Interface/GenEll/EllModuli.lean` の古い記述を訂正した）。

★`.src` は条つき——指標には数えない。
-/

namespace ABC3.Found.GenEll

open Matrix
open scoped MatrixGroups
open Matrix.SpecialLinearGroup

/-! ## ★★★★★★★★★★安定な直線が無ければ非上三角 -/

/-- ★★★★★★★★★★**安定な直線が無ければ非上三角な行列がある**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★対偶が自明である——**全部が上三角なら `⟨e₁⟩` が安定**だからである。
★★これが原文の『by the portion of Lemma 3.7 concerning `l`-cyclic subgroup schemes,
it follows that the image of Galois in `GL₂(𝔽_l)` contains at least one matrix which is
not upper triangular』の群論の側である
（`l`-巡回部分群スキーム ⟷ `E[l]` の中の `Gal`-安定な直線）。 -/
theorem exists_nonUpper_of_no_stable_line {l : ℕ} [Fact (Nat.Prime l)]
    (H : Set (GL (Fin 2) (ZMod l)))
    (hno : ¬ ∃ v : Fin 2 → ZMod l, v ≠ 0 ∧ ∀ M ∈ H, ∃ c : ZMod l,
        (M : Matrix (Fin 2) (Fin 2) (ZMod l)).mulVec v = c • v) :
    ∃ M ∈ H, (M : Matrix (Fin 2) (Fin 2) (ZMod l)) 1 0 ≠ 0 := by
  by_contra hc
  push_neg at hc
  apply hno
  refine ⟨Pi.single 0 1, ?_, ?_⟩
  · intro h
    have h0 := congrFun h 0
    simp at h0
  · intro M hM
    refine ⟨(M : Matrix (Fin 2) (Fin 2) (ZMod l)) 0 0, ?_⟩
    ext i
    fin_cases i <;> simp [Matrix.mulVec, hc M hM]

/-! ## ★★★★★★★★★★★★★★★原文 p.20 の最終段（群論の側） -/

/-- ★★★★★★★★★★★★★★★**原文 p.20 の最終段**（群論の側）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

`GL₂(ℤ_l)` の**閉**部分群 `J` について、mod `l` 像が

* `α = (1 1 / 0 1)` を含み、かつ
* **安定な直線を持たない**

なら、**`SL₂(ℤ_l) ⊆ J`**。

★`Lemma 3.1, (iv)`（`Sl2Padic.lean`）に `exists_nonUpper_of_no_stable_line` を差すだけである。
★★★これで `Theorem 3.8` に残るのは**幾何・数論の側だけ**になった:
局所理論の行列表示（`α` が入ること）・`l`-巡回 ⟷ 安定直線の対応・`Lemma 3.7`・`torsionExt`。 -/
theorem sl2_of_alpha_of_no_stable_line (l : ℕ) [Fact (Nat.Prime l)] (h5 : 5 ≤ l)
    (J : Subgroup (GL (Fin 2) ℤ_[l]))
    (hclosed : IsClosed (J : Set (GL (Fin 2) ℤ_[l])))
    (hα : (toGL (upper (1 : ZMod l)) : GL (Fin 2) (ZMod l)) ∈ J.map (glRedPadic l))
    (hno : ¬ ∃ v : Fin 2 → ZMod l, v ≠ 0 ∧
      ∀ M ∈ (↑(J.map (glRedPadic l)) : Set (GL (Fin 2) (ZMod l))), ∃ c : ZMod l,
        (M : Matrix (Fin 2) (Fin 2) (ZMod l)).mulVec v = c • v)
    (g : SL(2, ℤ_[l])) : (toGL g : GL (Fin 2) ℤ_[l]) ∈ J := by
  refine lemma_3_1_iv l h5 J hclosed hα ?_ g
  obtain ⟨M, hM, hMne⟩ := exists_nonUpper_of_no_stable_line _ hno
  exact ⟨M, hM, hMne⟩

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def exists_nonUpper_of_no_stable_line.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(安定な直線が無ければ非上三角な行列がある)",
    sectionId := "genell-thm-3-8" }

def sl2_of_alpha_of_no_stable_line.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(原文 p.20 の最終段——群論の側)",
    sectionId := "genell-thm-3-8" }

def sl2_of_alpha_of_no_stable_line.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_1_iv(Lemma 3.1, (iv)——SL₂(ℤ_l) の持ち上げ)"
      (.inProject "ABC3" "ABC3.Found.GenEll.lemma_3_1_iv") 4,
    .otherPaper "[GenEll]"
      ("Lemma 3.2 の直前の局所理論(l ∤ v(q) なら mod l 像が α を含む)" ++
       "——★Tate 一意化と Lemma 3.2, (i) は閉じている。**行列表示が残る**") 8,
    .folklore
      ("l-巡回部分群スキーム ⟷ E[l] の中の Gal-安定な直線の対応。**残る**") 5,
    .otherPaper "[GenEll]" "Lemma 3.7(局所高さと l-巡回部分群スキーム)——★Prop 3.4 待ち" 8,
    .implicitStep
      ("★★★★★★測定(2026-08-29): 原文 p.20 の最終段は 3 つでできており、" ++
       "★3 番(Lemma 3.1, (iv))は **Sl2Padic.lean に実装済み**、" ++
       "★★2 番の**群論の側**は本ファイルで取れた(対偶が自明——" ++
       "全部が上三角なら ⟨e₁⟩ が安定)。" ++
       "★★★残るのは幾何・数論の側だけである: " ++
       "局所理論の行列表示・l-巡回 ⟷ 安定直線の対応・Lemma 3.7・torsionExt") 7,
    .implicitStep
      ("★★Galois 表現そのものは Found/GaloisRep/GalRep.lean の galRep で**構成済み**である" ++
       "(2026-08-29 に Interface/GenEll/EllModuli.lean の古い記述を訂正した)。" ++
       "★したがって Theorem 3.8 に『Serre の開像定理』は要らない") 6 ]

end ABC3.Found.GenEll
