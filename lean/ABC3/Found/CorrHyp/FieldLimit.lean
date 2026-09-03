/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.RingTheory.Adjoin.FG

/-!
# [CorrHyp] `Lemma 4.1` へ向けた第一歩 —— `K` を有限生成 `k`-部分環の直極限として見る

`Lemma 4.1`(原文 p.9)の証明は「it only takes finitely many equations to
define a curve or a correspondence」という一文で、`X_K`・`Z_K`・その間の
correspondence が、ある有限生成 `k`-部分環 `R ⊆ K` 上ですでに定義できる
(spreading out)ことを使う。これを mathlib の
`AlgebraicGeometry/AffineTransitionLimit.lean`(スキームを余濾過的な
極限として扱う一般論、`ResearchPaper/corrhyp-goal.md` §3 参照)に載せる
ための、最初の(スキーム論に触れない)代数的な足場をここに置く。

★`FuchsianGroup`(`ℂ` 上の解析的モデル)とは別建てのトラック——`Lemma 4.1`
は一般の `k` を扱うため `ℍ`/`SL(2,ℝ)` のモデルは使えない。 -/

namespace ABC3.Found.CorrHyp

/-- `K` の元 `x, y` は、ある有限生成 `k`-部分環に同時に属する
(`Algebra.adjoin k {x, y}` を取ればよい)。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_fg_subalgebra_mem_two {k K : Type*} [CommRing k] [CommRing K] [Algebra k K]
    (x y : K) : ∃ R : Subalgebra k K, R.FG ∧ x ∈ R ∧ y ∈ R := by
  classical
  refine ⟨Algebra.adjoin k {x, y}, ?_, ?_, ?_⟩
  · have := Subalgebra.fg_adjoin_finset (R := k) ({x, y} : Finset K)
    simpa using this
  · exact Algebra.subset_adjoin (by simp)
  · exact Algebra.subset_adjoin (by simp)

/-- 有限生成 `k`-部分環は和(生成する部分環)を取っても有限生成
(`Subalgebra.FG.sup`、mathlib から無条件)。`Lemma 4.1` の証明が扱う
「有限個の方程式の係数をすべて含む1つの有限生成部分環」を作るのに使う
——複数の有限生成部分環を合わせても、なお有限生成のまま。 -/
theorem fg_sup {k K : Type*} [CommRing k] [CommRing K] [Algebra k K]
    (R S : Subalgebra k K) (hR : R.FG) (hS : S.FG) : (R ⊔ S).FG :=
  hR.sup hS

/-- `exists_fg_subalgebra_mem_two` の一般化: `K` の**任意有限個**の元は、
ある有限生成 `k`-部分環に同時に属する。原文の「it only takes finitely
many equations to define a curve or a correspondence」——`X_K`・`Z_K`・
その間の correspondence を定義する有限個の方程式の係数をまとめて
1つの有限生成部分環に落とす、という段の直接の道具になる。

★**sorry 無し**。標準3公理のみ。`exists_fg_subalgebra_mem_two` はこれの
特殊例なので、以後はこちらを使う。 -/
theorem exists_fg_subalgebra_mem_finset {k K : Type*} [CommRing k] [CommRing K] [Algebra k K]
    (s : Finset K) : ∃ R : Subalgebra k K, R.FG ∧ ∀ x ∈ s, x ∈ R := by
  refine ⟨Algebra.adjoin k (s : Set K), Subalgebra.fg_adjoin_finset s, ?_⟩
  intro x hx
  exact Algebra.subset_adjoin hx

/-- `K` の有限生成 `k`-部分環全体(包含で前順序)。`AffineTransitionLimit.lean`
の余濾過的な図式 `D : I ⥤ Scheme` の添字圏 `I` の候補——`Spec R`(`R` を
ここの元とする)たちが `Spec K` への遷移射を持つ図式になる。 -/
def FgSubalgebra (k K : Type*) [CommRing k] [CommRing K] [Algebra k K] : Type _ :=
  {R : Subalgebra k K // R.FG}

instance {k K : Type*} [CommRing k] [CommRing K] [Algebra k K] :
    Preorder (FgSubalgebra k K) :=
  Subtype.preorder _

/-- `FgSubalgebra k K` は有向集合をなす(`fg_sup` から——任意の2つの
有限生成部分環は、その和(なお有限生成)を共通の上界に持つ)。

★**sorry 無し**。標準3公理のみ。`AffineTransitionLimit.lean` の理論が
要求する「`I` が余濾過的(`IsCofiltered`)」の核となる性質——`IsDirected`
から `IsCofilteredOrEmpty`/`IsCofiltered`(前順序を圏と見たとき)が従う。 -/
instance {k K : Type*} [CommRing k] [CommRing K] [Algebra k K] :
    IsDirected (FgSubalgebra k K) (· ≤ ·) := by
  constructor
  intro ⟨R, hR⟩ ⟨S, hS⟩
  refine ⟨⟨R ⊔ S, fg_sup R S hR hS⟩, ?_, ?_⟩
  · show R ≤ R ⊔ S
    exact le_sup_left
  · show S ≤ R ⊔ S
    exact le_sup_right

end ABC3.Found.CorrHyp
