/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.RingTheory.Adjoin.FG
import Mathlib.Algebra.Category.Ring.Basic
import Mathlib.CategoryTheory.Category.Preorder
import Mathlib.AlgebraicGeometry.Scheme

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

/-- `K` はその有限生成 `k`-部分環全体の(有向)合併に等しい——
`K = ⋃ R`。双対を取れば `Spec K = lim Spec R`(`AffineTransitionLimit.lean`
が要求する形)になる、という主張の代数側の核。

★**sorry 無し**。標準3公理のみ。 -/
theorem iSup_fgSubalgebra_eq_top (k K : Type*) [CommRing k] [CommRing K] [Algebra k K] :
    (⨆ R : FgSubalgebra k K, R.1) = (⊤ : Subalgebra k K) := by
  apply le_antisymm le_top
  intro x _
  classical
  obtain ⟨R, hR, hx⟩ := exists_fg_subalgebra_mem_finset (k := k) ({x} : Finset K)
  exact le_iSup (fun R : FgSubalgebra k K => R.1) (⟨R, hR⟩ : FgSubalgebra k K)
    (hx x (Finset.mem_singleton_self x))

open CategoryTheory in
/-- `FgSubalgebra k K` を(前順序を薄い圏と見て)`CommRingCat` への図式にする
関手——`R ↦ R`(環として)、`R ≤ S` を包含環準同型 `R ↪ S` に送る。
`AffineTransitionLimit.lean` が要求する `D : I ⥤ Scheme` の手前、
`CommRingCat` レベルの図式(`Spec` と合成すれば `Scheme` の図式になる)。

★**sorry 無し**。標準3公理のみ。関手則(`map_id`/`map_comp`)は
`Subalgebra.inclusion` が包含写像そのものであることから `rfl` で落ちる。 -/
noncomputable def toRingCat (k K : Type*) [CommRing k] [CommRing K] [Algebra k K] :
    FgSubalgebra k K ⥤ CommRingCat where
  obj R := CommRingCat.of R.1
  map {R S} h := CommRingCat.ofHom (Subalgebra.inclusion (leOfHom h)).toRingHom
  map_id R := by
    apply CommRingCat.hom_ext
    ext x
    rfl
  map_comp {R S T} f g := by
    apply CommRingCat.hom_ext
    ext x
    rfl

open AlgebraicGeometry CategoryTheory in
/-- `toRingCat` を `Scheme.Spec`(`CommRingCatᵒᵖ ⥤ Scheme`)と合成した図式。
`FgSubalgebra k K` は(`⊔` で)**filtered**(`R ≤ S` の向きで合流点を持つ)
なので、その逆圏 `(FgSubalgebra k K)ᵒᵖ` が **cofiltered**——
`AffineTransitionLimit.lean` の `D : I ⥤ Scheme` が要求する形そのもの。
`R ↦ Spec R`、`R ≤ S` を `Spec S → Spec R` に送る。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def toSchemeDiagram (k K : Type*) [CommRing k] [CommRing K] [Algebra k K] :
    (FgSubalgebra k K)ᵒᵖ ⥤ Scheme :=
  (toRingCat k K).op ⋙ Scheme.Spec

/- ★★次の一手(未着手): `(FgSubalgebra k K)ᵒᵖ` が実際に
`CategoryTheory.IsCofiltered` のインスタンスを持つことを示す
(`IsDirected` から `IsFilteredOrEmpty`/`IsFiltered` を経由するはずだが、
`infer_instance` では自動的に繋がらないことを確認した——mathlib に
「有向前順序は filtered」という直接の instance が見当たらず、
`IsFiltered` の生の条件(cocone の存在・平行射の coequalize)から
手で組む必要がある、2026-09-04 実測)。これが済めば
`AffineTransitionLimit.lean` の spreading-out 定理群を `toSchemeDiagram`
に直接適用できる。 -/

end ABC3.Found.CorrHyp
