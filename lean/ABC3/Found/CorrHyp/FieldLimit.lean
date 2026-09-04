/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.RingTheory.Adjoin.FG
import Mathlib.Algebra.Category.Ring.Basic
import Mathlib.CategoryTheory.Category.Preorder
import Mathlib.AlgebraicGeometry.Scheme
import Mathlib.AlgebraicGeometry.GammaSpecAdjunction
import Mathlib.AlgebraicGeometry.AffineTransitionLimit
import Mathlib.Algebra.Polynomial.AlgebraMap
import Mathlib.Algebra.Polynomial.Monic
import Mathlib.RingTheory.Etale.StandardEtale
import Mathlib.RingTheory.TensorProduct.MvPolynomial
import Mathlib.RingTheory.TensorProduct.Quotient
import Mathlib.RingTheory.Unramified.LocalStructure

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

/-- **有限個の `FgSubalgebra k K` は共通の上界を持つ**——`IsDirected` を
`Finset` 上の有限帰納法で拡張しただけ(2つずつ合流させる)。`Lemma 4.1`
の「複数のアフィン片・複数の標準エタール片を横断した細分段階の合流」
(`corrhyp-goal.md` §4)で使う——有限個の局所降下先 `R_i` を1つの共通
段階へまとめる。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_fgSubalgebra_upperBound {k K : Type} [CommRing k] [CommRing K] [Algebra k K]
    {ι : Type} (t : Finset ι) (R : ι → FgSubalgebra k K) :
    ∃ R' : FgSubalgebra k K, ∀ i ∈ t, R i ≤ R' := by
  classical
  induction t using Finset.induction with
  | empty => exact ⟨⟨⊥, Subalgebra.fg_bot⟩, by simp⟩
  | insert a s ha ih =>
    obtain ⟨R', hR'⟩ := ih
    obtain ⟨R'', hR1, hR2⟩ := (instIsDirectedFgSubalgebraLe (k := k) (K := K)).directed (R a) R'
    refine ⟨R'', fun i hi => ?_⟩
    rcases Finset.mem_insert.mp hi with rfl | hi
    · exact hR1
    · exact le_trans (hR' i hi) hR2

/-- `exists_fgSubalgebra_upperBound`の2引数版(`Bool`添字の`Finset.univ`
経由)——「2つのRレベル段階の共通の上界」を毎回`Finset`/添字型を
書かずに直接呼べる、使い勝手のための特殊化。 -/
theorem exists_fgSubalgebra_upperBound2 {k K : Type} [CommRing k] [CommRing K] [Algebra k K]
    (R S : FgSubalgebra k K) : ∃ R', R ≤ R' ∧ S ≤ R' := by
  obtain ⟨R', hR'⟩ := exists_fgSubalgebra_upperBound (Finset.univ : Finset Bool)
    (fun b => if b then R else S)
  exact ⟨R', hR' true (Finset.mem_univ _), hR' false (Finset.mem_univ _)⟩

open CategoryTheory in
/-- `FgSubalgebra k K` は filtered(`IsDirected` から——2つの対象は `⊔` を
共通の余錐に持ち、薄い圏なので平行射の coequalize は自明)。

★**sorry 無し**。標準3公理のみ。 -/
instance instIsFilteredFgSubalgebra {k K : Type*} [CommRing k] [CommRing K] [Algebra k K] :
    IsFiltered (FgSubalgebra k K) where
  cocone_objs := fun ⟨R, hR⟩ ⟨S, hS⟩ =>
    ⟨⟨R ⊔ S, fg_sup R S hR hS⟩, homOfLE (le_sup_left (a := R) (b := S)),
      homOfLE (le_sup_right (a := R) (b := S)), trivial⟩
  cocone_maps {X Y} f g := ⟨Y, 𝟙 Y, Subsingleton.elim _ _⟩
  nonempty := ⟨⟨⊥, Subalgebra.fg_bot⟩⟩

open CategoryTheory in
/-- `(FgSubalgebra k K)ᵒᵖ` は cofiltered——`AffineTransitionLimit.lean` の
`D : I ⥤ Scheme` が要求する `I` の性質そのもの(mathlib の
`isCofilteredOrEmpty_op_of_isFilteredOrEmpty` から無条件)。

★**sorry 無し**。標準3公理のみ。これで前回のコミットで「次の一手」として
残していた欠落が埋まった——`toSchemeDiagram` に spreading-out 定理群を
直接適用できる状態になった。 -/
instance {k K : Type*} [CommRing k] [CommRing K] [Algebra k K] :
    IsCofiltered (FgSubalgebra k K)ᵒᵖ := by
  infer_instance

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

/-!
## `Spec K` が `toSchemeDiagram` の極限であること

前回の「次の一手」——`toSchemeDiagram` の極限錐を実際に構成し、頂点が `Spec K` に
一致することを示す——をここで埋める。道筋は「環側で余極限を作り、`op` を取って
`Scheme.Spec`(`Γ ⊣ Spec` の右随伴——極限を保つ)で像を送る」。

1. `toRingCat k K : FgSubalgebra k K ⥤ CommRingCat` の余錐 `toRingCatCocone`
   (頂点 `K`、各射は部分環の包含 `R ↪ K`)を作る。
2. その余錐が余極限であること(`isColimitToRingCatCocone`)を示す——
   これが本質的な内容: 余錐の値は代表元 `R` の取り方に依らない
   (`cocone_ι_congr'`、`FgSubalgebra k K` が有向集合をなすことから)ので、
   `K` の各元を「ある有限生成部分環に属する元」として持ち上げる環準同型
   `desc`(`toRingCatDesc`)が矛盾なく定義できる——mathlib に
   `Subalgebra.iSupLift`(`AlgHom` 版)はあるが `k`-線形性を要求しない
   `RingHom` 版が無かったので、`exists_fg_subalgebra_mem_finset` から
   手で組んだ。
3. `IsColimit.op` で余極限を `(toRingCatCocone k K).op : Cone (toRingCat k K).op`
   の極限に変換し、`Scheme.Spec` が極限を保つこと(`Γ ⊣ Spec` の右随伴、
   `ΓSpec.adjunction.rightAdjoint_preservesLimits`)で `Scheme.Spec` の像へ運ぶ。
   `(toRingCat k K).op ⋙ Scheme.Spec` は定義から `toSchemeDiagram k K` そのもの
   なので、`specKCone k K : Cone (toSchemeDiagram k K)` が得られ、その頂点は
   定義上 `Scheme.Spec.obj (op (CommRingCat.of K))` = `Spec K`。

★配管メモ: `FgSubalgebra k K := {R : Subalgebra k K // R.FG}` が `def`(非簡約)
であるため、`Category.comp_id`/`Functor.const_obj_map` 等を Preorder-as-category
の文脈で simp/rw しようとすると「not type-correct under the `instances`
transparency level」で止まる——mathlib の `AffineTransitionLimit.lean` 自身が
同種の場面で使っている `set_option backward.isDefEq.respectTransparency false`
が直接効く(`tools/lean-idioms.md` に追記)。
-/

open CategoryTheory Limits in
set_option backward.isDefEq.respectTransparency false in
/-- 余錐の値は「どの `R` を代表元に選ぶか」に依らない——`R ≤ T` のとき、余錐の
自然性から `T` 上で評価した値と `R` 上で評価した値が一致する。

★**sorry 無し**。標準3公理のみ。 -/
theorem cocone_ι_congr {k K : Type*} [CommRing k] [CommRing K] [Algebra k K]
    (s : Cocone (toRingCat k K)) {R T : FgSubalgebra k K} (h : R ≤ T) (x : K)
    (hR : x ∈ R.1) : s.ι.app T ⟨x, h hR⟩ = s.ι.app R ⟨x, hR⟩ := by
  have h1 := s.ι.naturality (homOfLE h)
  have h1' : (toRingCat k K).map (homOfLE h) ≫ s.ι.app T = s.ι.app R := by
    convert h1 using 2
    exact (Category.comp_id _).symm
  rw [← h1', CommRingCat.comp_apply]
  rfl

open CategoryTheory Limits in
/-- `cocone_ι_congr` の両側版: `R` と `S` が(有向性から)共通の上界 `T` を持つことを
使って、任意の2つの代表元での評価が一致することを示す。

★**sorry 無し**。標準3公理のみ。 -/
theorem cocone_ι_congr' {k K : Type*} [CommRing k] [CommRing K] [Algebra k K]
    (s : Cocone (toRingCat k K)) {R S : FgSubalgebra k K} (x : K) (hR : x ∈ R.1)
    (hS : x ∈ S.1) : s.ι.app R ⟨x, hR⟩ = s.ι.app S ⟨x, hS⟩ := by
  obtain ⟨T, hRT, hST⟩ := (instIsDirectedFgSubalgebraLe (k := k) (K := K)).directed R S
  rw [← cocone_ι_congr s hRT x hR, ← cocone_ι_congr s hST x hS]

open CategoryTheory Limits in
/-- `K` の元 `x` を、それを含むある有限生成部分環 `R`(`exists_fg_subalgebra_mem_finset`
で選ぶ)上での余錐の値として評価する——これが `desc` 環準同型の台になる。

★**sorry 無し**。標準3公理のみ。`noncomputable`(`Exists.choose` を使う)。 -/
noncomputable def descVal {k K : Type*} [CommRing k] [CommRing K] [Algebra k K]
    (s : Cocone (toRingCat k K)) (x : K) : s.pt :=
  s.ι.app ⟨(exists_fg_subalgebra_mem_finset (k := k) ({x} : Finset K)).choose,
      (exists_fg_subalgebra_mem_finset (k := k) ({x} : Finset K)).choose_spec.1⟩
    ⟨x, (exists_fg_subalgebra_mem_finset (k := k) ({x} : Finset K)).choose_spec.2
      x (Finset.mem_singleton_self x)⟩

open CategoryTheory Limits in
/-- `descVal` は実際にはどの `R`(`x` を含むもの)上で評価しても同じ値になる
(`cocone_ι_congr'` の言い換え)。以後、`descVal` の計算はすべてこれを経由する。

★**sorry 無し**。標準3公理のみ。 -/
theorem descVal_spec {k K : Type*} [CommRing k] [CommRing K] [Algebra k K]
    (s : Cocone (toRingCat k K)) {R : FgSubalgebra k K} (x : K) (hx : x ∈ R.1) :
    descVal s x = s.ι.app R ⟨x, hx⟩ := by
  unfold descVal
  apply cocone_ι_congr'

open CategoryTheory Limits in
/-- `descVal` を環準同型として束ねたもの——`K` の任意有限個の元は共通の有限生成
部分環に同時に属する(`exists_fg_subalgebra_mem_finset`)ので、`+`・`*` の
保存は「共通の代表元 `R` を選んで、`s.ι.app R` 自身が環準同型であることを使う」
だけで閉じる。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def toRingCatDesc {k K : Type*} [CommRing k] [CommRing K] [Algebra k K]
    (s : Cocone (toRingCat k K)) : K →+* s.pt where
  toFun := descVal s
  map_one' := by
    obtain ⟨R, hRfg, hRmem⟩ := exists_fg_subalgebra_mem_finset (k := k) (∅ : Finset K)
    have h1 : (1 : K) ∈ R := R.one_mem
    rw [descVal_spec s (R := ⟨R, hRfg⟩) 1 h1]
    show (ConcreteCategory.hom (s.ι.app ⟨R, hRfg⟩)) (1 : R) = 1
    exact map_one _
  map_zero' := by
    obtain ⟨R, hRfg, hRmem⟩ := exists_fg_subalgebra_mem_finset (k := k) (∅ : Finset K)
    have h1 : (0 : K) ∈ R := R.zero_mem
    rw [descVal_spec s (R := ⟨R, hRfg⟩) 0 h1]
    show (ConcreteCategory.hom (s.ι.app ⟨R, hRfg⟩)) (0 : R) = 0
    exact map_zero _
  map_mul' := fun x y => by
    classical
    obtain ⟨R, hRfg, hRmem⟩ := exists_fg_subalgebra_mem_finset (k := k) ({x, y} : Finset K)
    have hx : x ∈ R := hRmem x (by simp)
    have hy : y ∈ R := hRmem y (by simp)
    have hxy : x * y ∈ R := R.mul_mem hx hy
    rw [descVal_spec s (R := ⟨R, hRfg⟩) x hx, descVal_spec s (R := ⟨R, hRfg⟩) y hy,
      descVal_spec s (R := ⟨R, hRfg⟩) (x * y) hxy]
    show (ConcreteCategory.hom (s.ι.app ⟨R, hRfg⟩)) (⟨x, hx⟩ * ⟨y, hy⟩ : R) =
      (ConcreteCategory.hom (s.ι.app ⟨R, hRfg⟩)) ⟨x, hx⟩ *
        (ConcreteCategory.hom (s.ι.app ⟨R, hRfg⟩)) ⟨y, hy⟩
    exact map_mul _ _ _
  map_add' := fun x y => by
    classical
    obtain ⟨R, hRfg, hRmem⟩ := exists_fg_subalgebra_mem_finset (k := k) ({x, y} : Finset K)
    have hx : x ∈ R := hRmem x (by simp)
    have hy : y ∈ R := hRmem y (by simp)
    have hxy : x + y ∈ R := R.add_mem hx hy
    rw [descVal_spec s (R := ⟨R, hRfg⟩) x hx, descVal_spec s (R := ⟨R, hRfg⟩) y hy,
      descVal_spec s (R := ⟨R, hRfg⟩) (x + y) hxy]
    show (ConcreteCategory.hom (s.ι.app ⟨R, hRfg⟩)) (⟨x, hx⟩ + ⟨y, hy⟩ : R) =
      (ConcreteCategory.hom (s.ι.app ⟨R, hRfg⟩)) ⟨x, hx⟩ +
        (ConcreteCategory.hom (s.ι.app ⟨R, hRfg⟩)) ⟨y, hy⟩
    exact map_add _ _ _

open CategoryTheory Limits in
/-- `toRingCat k K` の標準余錐: 頂点 `K`、各射は部分環の包含 `R ↪ K`。 -/
noncomputable def toRingCatCocone (k K : Type*) [CommRing k] [CommRing K] [Algebra k K] :
    Cocone (toRingCat k K) where
  pt := CommRingCat.of K
  ι :=
  { app := fun R => CommRingCat.ofHom R.1.subtype
    naturality := by
      intro R S h
      apply CommRingCat.hom_ext
      ext x
      rfl }

open CategoryTheory Limits in
set_option backward.isDefEq.respectTransparency false in
private theorem toRingCatCocone_fac {k K : Type*} [CommRing k] [CommRing K] [Algebra k K]
    (s : Cocone (toRingCat k K)) (R : FgSubalgebra k K) :
    (toRingCatCocone k K).ι.app R ≫ CommRingCat.ofHom (toRingCatDesc s) = s.ι.app R := by
  apply CommRingCat.hom_ext
  ext x
  rw [CommRingCat.comp_apply]
  show toRingCatDesc s (x.1 : K) = s.ι.app R x
  show descVal s (x.1 : K) = s.ι.app R x
  exact descVal_spec s x.1 x.2

open CategoryTheory Limits in
set_option backward.isDefEq.respectTransparency false in
private theorem toRingCatCocone_uniq {k K : Type*} [CommRing k] [CommRing K] [Algebra k K]
    (s : Cocone (toRingCat k K)) (m : (toRingCatCocone k K).pt ⟶ s.pt)
    (hm : ∀ R : FgSubalgebra k K, (toRingCatCocone k K).ι.app R ≫ m = s.ι.app R) :
    m = CommRingCat.ofHom (toRingCatDesc s) := by
  apply CommRingCat.hom_ext
  ext x
  classical
  obtain ⟨R, hRfg, hRmem⟩ := exists_fg_subalgebra_mem_finset (k := k) ({x} : Finset K)
  have hx : x ∈ R := hRmem x (Finset.mem_singleton_self x)
  have h1 := hm ⟨R, hRfg⟩
  have h2 := congrArg (ConcreteCategory.hom ·) h1
  have h3 := DFunLike.congr_fun h2 (⟨x, hx⟩ : R)
  rw [CommRingCat.comp_apply] at h3
  show (ConcreteCategory.hom m) x = toRingCatDesc s x
  show (ConcreteCategory.hom m) x = descVal s x
  rw [descVal_spec s (R := ⟨R, hRfg⟩) x hx]
  exact h3

open CategoryTheory Limits in
/-- `toRingCatCocone` は実際に余極限——`K` がその有限生成部分環全体の
(有向な)合併であることの、圏論的な言い換え。

★**sorry 無し**。標準3公理のみ。これが本ファイルの核心の結果。 -/
noncomputable def isColimitToRingCatCocone (k K : Type*) [CommRing k] [CommRing K]
    [Algebra k K] : IsColimit (toRingCatCocone k K) :=
  IsColimit.mk (fun s => CommRingCat.ofHom (toRingCatDesc s))
    (fun s R => toRingCatCocone_fac s R) (fun s m hm => toRingCatCocone_uniq s m hm)

open CategoryTheory Limits AlgebraicGeometry in
/-- `Scheme.Spec` は極限を保つ(`Γ ⊣ Spec` の右随伴であることから、一般の
「右随伴は極限を保つ」定理 `Adjunction.rightAdjoint_preservesLimits` で無条件)。 -/
noncomputable instance preservesLimitsSpec : PreservesLimitsOfSize Scheme.Spec :=
  ΓSpec.adjunction.rightAdjoint_preservesLimits

open CategoryTheory Limits AlgebraicGeometry in
/-- `toSchemeDiagram k K` の極限錐——頂点は定義から `Spec K`
(`specKCone_pt` 参照)。`isColimitToRingCatCocone` を `.op` で反転し、
`Scheme.Spec` で極限として送ったもの。 -/
noncomputable def specKCone (k K : Type*) [CommRing k] [CommRing K] [Algebra k K] :
    Cone (toSchemeDiagram k K) :=
  Scheme.Spec.mapCone (toRingCatCocone k K).op

open CategoryTheory Limits AlgebraicGeometry in
/-- `specKCone` の頂点は(定義から)`Spec K` そのもの。

★**sorry 無し**。`rfl` で閉じる——`toSchemeDiagram`・`specKCone` の構成が
まさにそうなるように組んである。 -/
theorem specKCone_pt (k K : Type*) [CommRing k] [CommRing K] [Algebra k K] :
    (specKCone k K).pt = Scheme.Spec.obj (Opposite.op (CommRingCat.of K)) := rfl

open CategoryTheory Limits AlgebraicGeometry in
/-- **`Spec K` は `toSchemeDiagram k K` の極限である**——本ファイルの目標。
これで `c := specKCone k K`・`hc := isLimit_specKCone k K` として、
`AffineTransitionLimit.lean` の spreading-out 定理群
(`Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation`/
`Scheme.exists_hom_comp_eq_comp_of_locallyOfFiniteType`)を直接呼べる状態になった。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def isLimit_specKCone (k K : Type*) [CommRing k] [CommRing K] [Algebra k K] :
    IsLimit (specKCone k K) :=
  isLimitOfPreserves Scheme.Spec (IsColimit.op (isColimitToRingCatCocone k K))

/-!
## `toSchemeDiagram` の側条件——`AffineTransitionLimit.lean` の定理群がそのまま使える形

`Scheme.exists_hom_comp_eq_comp_of_locallyOfFiniteType`/
`Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation` はどちらも
`[IsCofiltered I]`(済み)に加え、`∀ {i j} (f : i ⟶ j), IsAffineHom (D.map f)`・
`∀ i, CompactSpace (D.obj i)`・`∀ i, QuasiSeparatedSpace (D.obj i)` を要求する。
`D := toSchemeDiagram k K` の場合、これらは`D.obj i`・`D.map f`が常に
`Spec (…)`の形であることから mathlib の一般的なインスタンス(アフィンスキーム
は常にコンパクト・準分離、環準同型の `Spec` は常にアフィン射)だけで
**無条件に**満たされる——`toSchemeDiagram`(`Functor.comp`)越しには
インスタンス探索が `Spec` の形を直接見ないので、`show` で `Spec` の形に
戻してから `infer_instance` に渡す3つの instance をここに置く。 -/

open AlgebraicGeometry CategoryTheory in
instance toSchemeDiagram_isAffineHom {k K : Type*} [CommRing k] [CommRing K] [Algebra k K]
    {i j : (FgSubalgebra k K)ᵒᵖ} (f : i ⟶ j) :
    IsAffineHom ((toSchemeDiagram k K).map f) := by
  show IsAffineHom (Scheme.Spec.map ((toRingCat k K).op.map f))
  infer_instance

open AlgebraicGeometry CategoryTheory in
instance toSchemeDiagram_compactSpace {k K : Type*} [CommRing k] [CommRing K] [Algebra k K]
    (i : (FgSubalgebra k K)ᵒᵖ) : CompactSpace ((toSchemeDiagram k K).obj i) := by
  show CompactSpace (Scheme.Spec.obj ((toRingCat k K).op.obj i))
  infer_instance

open AlgebraicGeometry CategoryTheory in
instance toSchemeDiagram_quasiSeparatedSpace {k K : Type*} [CommRing k] [CommRing K]
    [Algebra k K] (i : (FgSubalgebra k K)ᵒᵖ) :
    QuasiSeparatedSpace ((toSchemeDiagram k K).obj i) := by
  show QuasiSeparatedSpace (Scheme.Spec.obj ((toRingCat k K).op.obj i))
  infer_instance

/-!
## 多項式の係数を有限段階へ降ろす

`Lemma 4.1` の構成的な降下(`corrhyp-goal.md` §4 に記録した見通し)の最初の
歯車: 有限エタール射の標準的な表示(`Algebra.StandardEtalePresentation`、
`f`・`g : Polynomial K` という2つの多項式データ)は、係数が有限個しかないので、
ある有限生成 `k`-部分環 `R` 上の多項式に必ず降ろせる——`monic` 性も
`algebraMap R K` の単射性から遺伝する。 -/

open Polynomial in
/-- `K` の多項式 `p` は、ある有限生成 `k`-部分環 `R` 上の多項式の像として書ける
——`p` の(有限個の)係数がすべてある `R` に同時に属することから。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_fg_subalgebra_polynomial {k K : Type*} [CommRing k] [CommRing K] [Algebra k K]
    (p : Polynomial K) :
    ∃ (R : FgSubalgebra k K) (p₀ : Polynomial R.1), p₀.map (algebraMap R.1 K) = p := by
  classical
  obtain ⟨R, hRfg, hRmem⟩ :=
    exists_fg_subalgebra_mem_finset (k := k) (p.support.image p.coeff)
  have hmem : ∀ x ∈ p.support, p.coeff x ∈ R :=
    fun x hx => hRmem (p.coeff x) (Finset.mem_image_of_mem p.coeff hx)
  refine ⟨⟨R, hRfg⟩, p.support.attach.sum (fun i =>
    Polynomial.monomial (i : ℕ) (⟨p.coeff (i : ℕ), hmem i i.2⟩ : R)), ?_⟩
  rw [Polynomial.map_sum]
  have key : ∀ i : {x // x ∈ p.support}, (Polynomial.monomial (i : ℕ)
      (⟨p.coeff (i : ℕ), hmem i i.2⟩ : R)).map (algebraMap ↥R K) =
      Polynomial.monomial (i : ℕ) (p.coeff (i : ℕ)) := by
    intro i
    rw [Polynomial.map_monomial]
    rfl
  simp_rw [key]
  rw [Finset.sum_attach p.support (fun i => Polynomial.monomial i (p.coeff i))]
  exact Polynomial.sum_monomial_eq p

open Polynomial in
/-- `exists_fg_subalgebra_polynomial` の2変数版——`p`・`q` の係数をまとめて
1つの有限生成部分環に降ろす(標準エタール表示の `f`・`g` を同時に降ろすのに使う)。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_fg_subalgebra_polynomial_pair {k K : Type*} [CommRing k] [CommRing K] [Algebra k K]
    (p q : Polynomial K) :
    ∃ (R : FgSubalgebra k K) (p₀ q₀ : Polynomial R.1),
      p₀.map (algebraMap R.1 K) = p ∧ q₀.map (algebraMap R.1 K) = q := by
  classical
  obtain ⟨R, hRfg, hRmem⟩ :=
    exists_fg_subalgebra_mem_finset (k := k) (p.support.image p.coeff ∪ q.support.image q.coeff)
  have hpmem : ∀ x ∈ p.support, p.coeff x ∈ R := fun x hx =>
    hRmem (p.coeff x) (Finset.mem_union_left _ (Finset.mem_image_of_mem p.coeff hx))
  have hqmem : ∀ x ∈ q.support, q.coeff x ∈ R := fun x hx =>
    hRmem (q.coeff x) (Finset.mem_union_right _ (Finset.mem_image_of_mem q.coeff hx))
  refine ⟨⟨R, hRfg⟩,
    p.support.attach.sum (fun i => Polynomial.monomial (i : ℕ) (⟨p.coeff (i : ℕ), hpmem i i.2⟩ : R)),
    q.support.attach.sum (fun i => Polynomial.monomial (i : ℕ) (⟨q.coeff (i : ℕ), hqmem i i.2⟩ : R)),
    ?_, ?_⟩
  · rw [Polynomial.map_sum]
    have key : ∀ i : {x // x ∈ p.support}, (Polynomial.monomial (i : ℕ)
        (⟨p.coeff (i : ℕ), hpmem i i.2⟩ : R)).map (algebraMap ↥R K) =
        Polynomial.monomial (i : ℕ) (p.coeff (i : ℕ)) := by
      intro i; rw [Polynomial.map_monomial]; rfl
    simp_rw [key]
    rw [Finset.sum_attach p.support (fun i => Polynomial.monomial i (p.coeff i))]
    exact Polynomial.sum_monomial_eq p
  · rw [Polynomial.map_sum]
    have key : ∀ i : {x // x ∈ q.support}, (Polynomial.monomial (i : ℕ)
        (⟨q.coeff (i : ℕ), hqmem i i.2⟩ : R)).map (algebraMap ↥R K) =
        Polynomial.monomial (i : ℕ) (q.coeff (i : ℕ)) := by
      intro i; rw [Polynomial.map_monomial]; rfl
    simp_rw [key]
    rw [Finset.sum_attach q.support (fun i => Polynomial.monomial i (q.coeff i))]
    exact Polynomial.sum_monomial_eq q

open Polynomial in
/-- `exists_fg_subalgebra_polynomial_pair` に「`p` が monic なら `p₀` も monic」を
足した版——`algebraMap R K`(部分環の包含)が単射であることから
(`Function.Injective.monic_map_iff`)。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_fg_subalgebra_polynomial_pair_monic {k K : Type*} [CommRing k] [CommRing K]
    [Algebra k K] (p q : Polynomial K) (hp : p.Monic) :
    ∃ (R : FgSubalgebra k K) (p₀ q₀ : Polynomial R.1),
      p₀.map (algebraMap R.1 K) = p ∧ q₀.map (algebraMap R.1 K) = q ∧ p₀.Monic := by
  obtain ⟨R, p₀, q₀, hp₀, hq₀⟩ := exists_fg_subalgebra_polynomial_pair (k := k) p q
  refine ⟨R, p₀, q₀, hp₀, hq₀, ?_⟩
  have hinj : Function.Injective (algebraMap R.1 K) := Subtype.coe_injective
  refine (hinj.monic_map_iff (p := p₀)).mpr ?_
  rw [hp₀]; exact hp

open Polynomial in
/-- `exists_fg_subalgebra_polynomial` の有限族版——有限個の多項式をまとめて
1つの有限生成部分環に同時に降ろす(`StandardEtalePair` の `cond` が要求する
4つの多項式 `f, g, p₁, p₂` を同時に降ろすのに使う)。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_fg_subalgebra_polynomial_family {k K : Type*} [CommRing k] [CommRing K]
    [Algebra k K] {ι : Type*} [Fintype ι] [DecidableEq K] (p : ι → Polynomial K) :
    ∃ (R : FgSubalgebra k K) (p₀ : ι → Polynomial R.1),
      ∀ i, (p₀ i).map (algebraMap R.1 K) = p i := by
  classical
  obtain ⟨R, hRfg, hRmem⟩ := exists_fg_subalgebra_mem_finset (k := k)
    (Finset.univ.biUnion (fun i => (p i).support.image (p i).coeff))
  have hmem : ∀ i (x : ℕ), x ∈ (p i).support → (p i).coeff x ∈ R := fun i x hx =>
    hRmem ((p i).coeff x) (Finset.mem_biUnion.mpr ⟨i, Finset.mem_univ i,
      Finset.mem_image_of_mem (p i).coeff hx⟩)
  refine ⟨⟨R, hRfg⟩, fun i => (p i).support.attach.sum (fun j =>
    Polynomial.monomial (j : ℕ) (⟨(p i).coeff (j : ℕ), hmem i j j.2⟩ : R)), fun i => ?_⟩
  rw [Polynomial.map_sum]
  have key : ∀ j : {x // x ∈ (p i).support}, (Polynomial.monomial (j : ℕ)
      (⟨(p i).coeff (j : ℕ), hmem i j j.2⟩ : R)).map (algebraMap ↥R K) =
      Polynomial.monomial (j : ℕ) ((p i).coeff (j : ℕ)) := by
    intro j; rw [Polynomial.map_monomial]; rfl
  simp_rw [key]
  rw [Finset.sum_attach (p i).support (fun j => Polynomial.monomial j ((p i).coeff j))]
  exact Polynomial.sum_monomial_eq (p i)

open Polynomial in
/-- `Algebra.StandardEtalePair` の `cond`(`f' * p₁ + f * p₂ = g^n` を満たす
`p₁, p₂, n` の存在)は、`f`・`g`・`p₁`・`p₂` の係数が有限個であることから、
ある有限生成 `k`-部分環上ですでに成り立つ形に降ろせる——4つの多項式を
`exists_fg_subalgebra_polynomial_family` で同時に降ろし、`Polynomial.map`
の単射性(`algebraMap R K` が単射なので)で等式そのものを引き戻す。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_fg_subalgebra_standardEtaleCond {k K : Type*} [CommRing k] [CommRing K]
    [Algebra k K] (f g p₁ p₂ : Polynomial K) (n : ℕ)
    (hcond : Polynomial.derivative f * p₁ + f * p₂ = g ^ n) :
    ∃ (R : FgSubalgebra k K) (f₀ g₀ p₁₀ p₂₀ : Polynomial R.1),
      f₀.map (algebraMap R.1 K) = f ∧ g₀.map (algebraMap R.1 K) = g ∧
      p₁₀.map (algebraMap R.1 K) = p₁ ∧ p₂₀.map (algebraMap R.1 K) = p₂ ∧
      Polynomial.derivative f₀ * p₁₀ + f₀ * p₂₀ = g₀ ^ n := by
  classical
  obtain ⟨R, q, hq⟩ := exists_fg_subalgebra_polynomial_family (k := k)
    (![f, g, p₁, p₂] : Fin 4 → Polynomial K)
  have hf : (q 0).map (algebraMap R.1 K) = f := by simpa using hq 0
  have hg : (q 1).map (algebraMap R.1 K) = g := by simpa using hq 1
  have hp1 : (q 2).map (algebraMap R.1 K) = p₁ := by simpa using hq 2
  have hp2 : (q 3).map (algebraMap R.1 K) = p₂ := by simpa using hq 3
  refine ⟨R, q 0, q 1, q 2, q 3, hf, hg, hp1, hp2, ?_⟩
  have hinj : Function.Injective (Polynomial.map (algebraMap R.1 K)) :=
    Polynomial.map_injective _ Subtype.coe_injective
  apply hinj
  rw [Polynomial.map_add, Polynomial.map_mul, Polynomial.map_mul, Polynomial.map_pow,
    ← Polynomial.derivative_map, hf, hg, hp1, hp2, hcond]

/-- **`Algebra.StandardEtalePair` は有限段階へ降ろせる**——`f`(monic)・`g`・
`cond` の証拠 `p₁, p₂, n` の係数がすべて有限個であることから、ある有限生成
`k`-部分環 `R` 上の `StandardEtalePair` の像として書ける。`Lemma 4.1` の
構成的降下(`corrhyp-goal.md` §4)の核となる部品——有限エタール射の
標準的表示が丸ごと有限段階に降りることを保証する。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_fg_subalgebra_standardEtalePair {k K : Type*} [CommRing k] [CommRing K]
    [Algebra k K] (P : StandardEtalePair K) :
    ∃ (R : FgSubalgebra k K) (P₀ : StandardEtalePair R.1),
      P₀.f.map (algebraMap R.1 K) = P.f ∧ P₀.g.map (algebraMap R.1 K) = P.g := by
  obtain ⟨p₁, p₂, n, hcond⟩ := P.cond
  obtain ⟨R, f₀, g₀, p₁₀, p₂₀, hf, hg, hp1, hp2, hcond₀⟩ :=
    exists_fg_subalgebra_standardEtaleCond (k := k) P.f P.g p₁ p₂ n hcond
  have hinj : Function.Injective (algebraMap R.1 K) := Subtype.coe_injective
  have hmonic : f₀.Monic := by
    refine (hinj.monic_map_iff (p := f₀)).mpr ?_
    rw [hf]; exact P.monic_f
  exact ⟨R, ⟨f₀, hmonic, g₀, ⟨p₁₀, p₂₀, n, hcond₀⟩⟩, hf, hg⟩

/-- `exists_fg_subalgebra_standardEtalePair` を mathlib の
`StandardEtalePair.map`(係数写像による像、`(P.map φ).f = P.f.map φ` 等)の
言葉で言い換えたもの——`P` は、ある有限生成部分環上の `P₀` の
`.map (algebraMap R K)` に**文字通り一致する**(`f`・`g` が一致すれば、
`monic_f`・`cond` は `Prop` なので自動的に一致する)。`P₀.Ring`
(`Polynomial (Polynomial R) ⧸ span {C f, X*C g - 1}`、mathlib の明示的な
構成)を使って `Z` の局所片を直接定義する際の橋渡しに使う。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_fg_subalgebra_standardEtalePair_map {k K : Type*} [CommRing k] [CommRing K]
    [Algebra k K] (P : StandardEtalePair K) :
    ∃ (R : FgSubalgebra k K) (P₀ : StandardEtalePair R.1), P₀.map (algebraMap R.1 K) = P := by
  obtain ⟨R, P₀, hf, hg⟩ := exists_fg_subalgebra_standardEtalePair (k := k) P
  refine ⟨R, P₀, ?_⟩
  obtain ⟨f, hfmonic, g, hcond⟩ := P
  obtain ⟨f₀, hfmonic₀, g₀, hcond₀⟩ := P₀
  simp only [StandardEtalePair.map] at hf hg ⊢
  subst hf hg
  rfl

/-!
## `StandardEtalePair.Ring` の base change

`P₀.Ring`(`Polynomial (Polynomial R) ⧸ span {C f, X*C g - 1}`、mathlib の
明示的な構成)が `K` へ係数拡大すると `(P₀.map φ).Ring` に一致することを示す
——`StandardEtalePair.equivMvPolynomialQuotient`(`P.Ring` を `MvPolynomial
(Fin 2) R` の商として書き直す)・`MvPolynomial.algebraTensorAlgEquiv`
(多変数多項式環の base change)・`Algebra.TensorProduct.tensorQuotientEquiv`
(商環の base change)を合成する。核心は「生成元の像が対応する生成元に写る」
ことの確認(`Bivariate_equivMvPolynomial_map`)。 -/

open Polynomial MvPolynomial Algebra.TensorProduct in
/-- `Polynomial.Bivariate.equivMvPolynomial` は係数の写像と可換
(`R[X][X] → MvPolynomial (Fin 2) R` という同型が、係数環の写像 `φ` と
両立する)——両辺とも `R[X][X] → MvPolynomial (Fin 2) S` への環準同型なので
`Polynomial.ringHom_ext'` を2回(外側の変数・内側の変数)適用して確認する。

★**sorry 無し**。標準3公理のみ。 -/
theorem Bivariate_equivMvPolynomial_map {R S : Type*} [CommRing R] [CommRing S] (φ : R →+* S)
    (p : R[X][X]) :
    (Bivariate.equivMvPolynomial S) (p.map (Polynomial.mapRingHom φ)) =
      MvPolynomial.map φ (Bivariate.equivMvPolynomial R p) := by
  have hF : (RingHom.comp (Bivariate.equivMvPolynomial S).toRingHom
        (Polynomial.mapRingHom (Polynomial.mapRingHom φ))) =
      (RingHom.comp (MvPolynomial.map φ) (Bivariate.equivMvPolynomial R).toRingHom) := by
    apply Polynomial.ringHom_ext'
    · apply Polynomial.ringHom_ext'
      · ext r
        simp [Bivariate.equivMvPolynomial_C_C]
      · simp [Bivariate.equivMvPolynomial_C_X]
    · simp [Bivariate.equivMvPolynomial_X]
  exact RingHom.congr_fun hF p

open Polynomial MvPolynomial Algebra.TensorProduct in
/-- `P₀.equivMvPolynomialQuotient` が使うイデアル(`{C f, X*C g - 1}` の
span)は、`includeRight` で `K ⊗[R] MvPolynomial (Fin 2) R` へ送ってから
`algebraTensorAlgEquiv R K`(`≃ₐ[K] MvPolynomial (Fin 2) K`)で運ぶと、
`(P₀.map φ)` の対応するイデアルに**文字通り一致する**——`Ideal.map_span` で
生成元の像を計算し、`Bivariate_equivMvPolynomial_map` で係数写像との可換性を
使う。

★**sorry 無し**。標準3公理のみ。 -/
theorem standardEtalePair_ring_baseChange {R K : Type*} [CommRing R] [CommRing K] [Algebra R K]
    (P₀ : StandardEtalePair R) :
    Ideal.map (algebraTensorAlgEquiv R K).toAlgHom.toRingHom
        (Ideal.map (includeRight (R := R) (A := K) (B := MvPolynomial (Fin 2) R))
          (Ideal.span {Bivariate.equivMvPolynomial R (Polynomial.C P₀.f),
            Bivariate.equivMvPolynomial R (Polynomial.X * Polynomial.C P₀.g - 1)})) =
      Ideal.span {Bivariate.equivMvPolynomial K (Polynomial.C (P₀.map (algebraMap R K)).f),
        Bivariate.equivMvPolynomial K
          (Polynomial.X * Polynomial.C (P₀.map (algebraMap R K)).g - 1)} := by
  rw [show (Ideal.map (includeRight (R := R) (A := K) (B := MvPolynomial (Fin 2) R))
          (Ideal.span {Bivariate.equivMvPolynomial R (Polynomial.C P₀.f),
            Bivariate.equivMvPolynomial R (Polynomial.X * Polynomial.C P₀.g - 1)}) :
          Ideal (TensorProduct R K (MvPolynomial (Fin 2) R))) =
      Ideal.map (includeRight (R := R) (A := K) (B := MvPolynomial (Fin 2) R)).toRingHom
          (Ideal.span {Bivariate.equivMvPolynomial R (Polynomial.C P₀.f),
            Bivariate.equivMvPolynomial R (Polynomial.X * Polynomial.C P₀.g - 1)}) from rfl]
  rw [Ideal.map_map, Ideal.map_span]
  congr 1
  simp only [Set.image_insert_eq, Set.image_singleton, StandardEtalePair.map_f,
    StandardEtalePair.map_g]
  congr 2
  · show (algebraTensorAlgEquiv R K)
      (includeRight (Bivariate.equivMvPolynomial R (Polynomial.C P₀.f))) = _
    rw [show includeRight (Bivariate.equivMvPolynomial R (Polynomial.C P₀.f)) =
      (1 : K) ⊗ₜ[R] (Bivariate.equivMvPolynomial R (Polynomial.C P₀.f)) from rfl,
      algebraTensorAlgEquiv_tmul, one_smul,
      ← Bivariate_equivMvPolynomial_map]
    congr 1
    simp [Polynomial.map_C]
  · show (algebraTensorAlgEquiv R K) (includeRight (Bivariate.equivMvPolynomial R
        (Polynomial.X * Polynomial.C P₀.g - 1))) = _
    rw [show includeRight (Bivariate.equivMvPolynomial R
        (Polynomial.X * Polynomial.C P₀.g - 1)) =
      (1 : K) ⊗ₜ[R] (Bivariate.equivMvPolynomial R (Polynomial.X * Polynomial.C P₀.g - 1))
        from rfl,
      algebraTensorAlgEquiv_tmul, one_smul,
      ← Bivariate_equivMvPolynomial_map]
    congr 1
    simp [Polynomial.map_sub, Polynomial.map_mul, Polynomial.map_C, Polynomial.map_X]

open Polynomial MvPolynomial Algebra.TensorProduct in
/-- **`StandardEtalePair.Ring` は base change と可換**——`P₀.Ring` を `K` へ
係数拡大すると `(P₀.map (algebraMap R K)).Ring` に(`K`-代数として)同型。
`Lemma 4.1` の構成的降下(`corrhyp-goal.md` §4)の核心部品:
`exists_fg_subalgebra_standardEtalePair_map` で得た有限段階の `P₀` から
`Z` の局所片 `Spec P₀.Ring` を作れば、この同型がまさに「その base change が
元の有限エタール多元環に一致する」ことを保証する。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def standardEtalePairRingBaseChange {R K : Type*} [CommRing R] [CommRing K]
    [Algebra R K] (P₀ : StandardEtalePair R) :
    TensorProduct R K P₀.Ring ≃ₐ[K] (P₀.map (algebraMap R K)).Ring :=
  (Algebra.TensorProduct.congr AlgEquiv.refl P₀.equivMvPolynomialQuotient).trans <|
    (tensorQuotientEquiv K (MvPolynomial (Fin 2) R) K
      (Ideal.span {Bivariate.equivMvPolynomial R (Polynomial.C P₀.f),
        Bivariate.equivMvPolynomial R (Polynomial.X * Polynomial.C P₀.g - 1)})).trans <|
    (Ideal.quotientEquivAlg _ _ (algebraTensorAlgEquiv R K)
      (standardEtalePair_ring_baseChange P₀).symm).trans
    (P₀.map (algebraMap R K)).equivMvPolynomialQuotient.symm

/- ★★次の一手(未着手): `Lemma 4.1` 本体へ——`standardEtalePairRingBaseChange`
により「有限段階の標準エタール表示 `P₀` から作った `Z` の局所片
`Spec P₀.Ring` の base change が、元の `K`-代数に(標準的な同型で)一致する」
という `Lemma 4.1` の構成的降下の核心部品が完成した。残るのは
(a) 一般の有限エタール多元環は「至る所で」標準エタールとは限らないため、
`Z_K` 全体のアフィン開被覆(`Scheme.exists_isOpenCover_and_isAffine`)の
**各片**をさらに étale-locus のレベルで細分してからこの降下を適用すること、
(b) 各片の降下結果を
`Scheme.exists_hom_hom_comp_eq_comp_of_locallyOfFiniteType`(遷移射の一致の
降下)で貼り合わせて `Z` 全体を構成すること。`corrhyp-goal.md` §4 に記録した
組み立て方の続き。 -/

/-!
## étale 多元環は有限個の基本開集合上で標準エタールになる

`exists_extDiagram_finite_affine_descent`(`ExtLimit.lean`)で `Z_K` の
アフィン開被覆を有限段階へ降ろした後に要る次の一手——各アフィン片の上の
有限エタール多元環は「至る所で」標準エタールとは限らないため、
`Algebra.IsEtaleAt.exists_isStandardEtale`(mathlib既存、各点で局所的に
標準エタールになる)を各素点に適用し、`Spec S` の準コンパクト性で
有限個の基本開集合に絞り込む——`ExtLimit.lean` の
`Scheme.exists_finite_affineOpenCover` と同じ「点ごとの局所的な性質→
準コンパクト性で有限化」というパターン。 -/

/-- **étale な多元環は、有限個の基本開集合の上で標準エタールになる**——
各素点で `Algebra.IsEtaleAt.exists_isStandardEtale` を適用し、得られる
基本開被覆を `PrimeSpectrum S` の準コンパクト性(`PrimeSpectrum.
compactSpace`)で有限部分被覆に絞り込む。CorrHyp に依存しない一般的な事実。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_finite_standardEtaleCover (R S : Type) [CommRing R] [CommRing S] [Algebra R S]
    [Algebra.Etale R S] :
    ∃ (ι : Type) (t : Finset ι) (f : ι → S),
      (⨆ i ∈ t, PrimeSpectrum.basicOpen (f i)) = ⊤ ∧
      ∀ i ∈ t, Algebra.IsStandardEtale R (Localization.Away (f i)) := by
  have hfp : Algebra.FinitePresentation R S := inferInstance
  choose f hf hstd using
    fun Q : PrimeSpectrum S => Algebra.IsEtaleAt.exists_isStandardEtale (R := R) Q.asIdeal
  have hcov : ⋃ Q : PrimeSpectrum S,
      (PrimeSpectrum.basicOpen (f Q) : Set (PrimeSpectrum S)) = Set.univ := by
    apply Set.eq_univ_of_forall
    intro Q
    exact Set.mem_iUnion.mpr ⟨Q, hf Q⟩
  obtain ⟨t, ht⟩ := (isCompact_univ (X := PrimeSpectrum S)).elim_finite_subcover
    (fun Q => (PrimeSpectrum.basicOpen (f Q) : Set (PrimeSpectrum S)))
    (fun Q => (PrimeSpectrum.basicOpen (f Q)).isOpen) (by rw [hcov])
  refine ⟨PrimeSpectrum S, t, f, ?_, fun i _ => hstd i⟩
  apply TopologicalSpace.Opens.ext
  rw [TopologicalSpace.Opens.coe_iSup]
  simp only [TopologicalSpace.Opens.coe_top]
  apply Set.eq_univ_of_univ_subset
  intro x _
  have hx := ht (Set.mem_univ x)
  simpa using hx

/-! ## `StandardEtalePair` の `A ⊗[ℚ] ℝ` からの降下——generic flatness 回避版

`Lemma 4.1` の「1アフィン片の降下」(`corrhyp-goal.md` §4、`piecePullbackIso`
と組み合わせて使う)向けの専用版。上の `exists_fg_subalgebra_standardEtalePair`
(一般の `k K`)と同型の構成だが、`k := ℚ`(体)に固定し、`K` の代わりに
`A ⊗[ℚ] ℝ`(`A` は任意の `ℚ`-代数、`Ext X` のアフィン片の切断環)を使う。
`ℚ` は体なので `A ⊗[ℚ] -` が自動的に平坦(`Module.Flat` が `infer_instance`
一発)——generic flatness(EGA IV、mathlib に無い)が不要になる鍵。 -/

open scoped TensorProduct in
/-- `A ⊗[ℚ] ℝ` の有限個の元は、ある `R : FgSubalgebra ℚ ℝ` を使って
`A ⊗[ℚ] R.1` からの像として同時に書ける——`exists_fg_subalgebra_mem_finset`
を各元のテンソル分解(`TensorProduct.exists_finset`)の係数へ適用する。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_fg_subalgebra_tensor_finset (A : Type) [CommRing A] [Algebra ℚ A]
    (s : Finset (A ⊗[ℚ] ℝ)) :
    ∃ R : FgSubalgebra ℚ ℝ, ∀ x ∈ s, ∃ y : A ⊗[ℚ] R.1,
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)) y = x := by
  classical
  choose S hS using fun x : s => TensorProduct.exists_finset (R := ℚ) (M := A) (N := ℝ) x.1
  set allR : Finset ℝ := s.attach.biUnion (fun x => (S x).image Prod.snd) with hallR
  obtain ⟨R, hRfg, hRmem⟩ := exists_fg_subalgebra_mem_finset (k := ℚ) allR
  refine ⟨⟨R, hRfg⟩, ?_⟩
  intro x hx
  set xs : s := ⟨x, hx⟩ with hxs
  have hmem : ∀ i ∈ S xs, i.2 ∈ R := by
    intro i hi
    apply hRmem
    rw [hallR]
    exact Finset.mem_biUnion.mpr ⟨xs, s.mem_attach xs, Finset.mem_image_of_mem Prod.snd hi⟩
  refine ⟨∑ i ∈ (S xs).attach, i.1.1 ⊗ₜ[ℚ] (⟨i.1.2, hmem i.1 i.2⟩ : R), ?_⟩
  rw [map_sum]
  simp only [Algebra.TensorProduct.map_tmul, AlgHom.id_apply, Subalgebra.coe_val]
  rw [Finset.sum_attach (S xs) (fun i => i.1 ⊗ₜ[ℚ] i.2)]
  exact (hS xs).symm

open Polynomial in
open scoped TensorProduct in
/-- `exists_fg_subalgebra_tensor_finset` の有限族版(多項式の係数)。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_fg_subalgebra_tensor_polynomial_family (A : Type) [CommRing A] [Algebra ℚ A]
    {ι : Type} [Fintype ι] [DecidableEq (A ⊗[ℚ] ℝ)] (p : ι → Polynomial (A ⊗[ℚ] ℝ)) :
    ∃ (R : FgSubalgebra ℚ ℝ) (p₀ : ι → Polynomial (A ⊗[ℚ] R.1)),
      ∀ i, (p₀ i).map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom
        = p i := by
  classical
  obtain ⟨R, hR⟩ := exists_fg_subalgebra_tensor_finset A
    (Finset.univ.biUnion (fun i => (p i).support.image (p i).coeff))
  choose c hc using fun q : Σ i, (p i).support => hR ((p q.1).coeff q.2)
    (Finset.mem_biUnion.mpr ⟨q.1, Finset.mem_univ q.1, Finset.mem_image_of_mem (p q.1).coeff q.2.2⟩)
  refine ⟨R, fun i => (p i).support.attach.sum (fun j => Polynomial.monomial (j : ℕ) (c ⟨i, j⟩)),
    fun i => ?_⟩
  rw [Polynomial.map_sum]
  have key : ∀ j : {x // x ∈ (p i).support}, (Polynomial.monomial (j : ℕ) (c ⟨i, j⟩)).map
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom =
      Polynomial.monomial (j : ℕ) ((p i).coeff (j : ℕ)) := by
    intro j
    rw [Polynomial.map_monomial]
    exact congrArg (Polynomial.monomial (j : ℕ)) (hc ⟨i, j⟩)
  simp_rw [key]
  rw [Finset.sum_attach (p i).support (fun j => Polynomial.monomial j ((p i).coeff j))]
  exact Polynomial.sum_monomial_eq (p i)

open scoped TensorProduct in
/-- `A ⊗[ℚ] R.1 → A ⊗[ℚ] ℝ`(`R.1 ↪ ℝ` に沿った base change)は常に単射
——`ℚ` は体なので `Module.Flat ℚ A` が任意の `A` で自動的に成り立つ
(`infer_instance`)ことから、`Module.Flat.iff_lTensor_preserves_
injective_linearMap` で移す。generic flatness を要らなくする核心。

★**sorry 無し**。標準3公理のみ。 -/
theorem tensor_map_injective_of_flat (A : Type) [CommRing A] [Algebra ℚ A]
    (R : FgSubalgebra ℚ ℝ) :
    Function.Injective (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)) := by
  haveI : Module.Flat ℚ A := inferInstance
  have hinj : Function.Injective (Subalgebra.val R.1) := Subtype.coe_injective
  have hrT := Module.Flat.iff_lTensor_preserves_injective_linearMap.mp this
    (Subalgebra.val R.1).toLinearMap hinj
  have hcomm : ∀ z : A ⊗[ℚ] R.1, (Subalgebra.val R.1).toLinearMap.lTensor A z =
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)) z := by
    intro z
    induction z using TensorProduct.induction_on with
    | zero => simp
    | tmul a r => simp
    | add z1 z2 h1 h2 => simp [LinearMap.map_add, h1, h2]
  intro x y hxy
  apply hrT
  rw [hcomm, hcomm, hxy]

open Polynomial in
open scoped TensorProduct in
/-- `exists_fg_subalgebra_standardEtaleCond` のtensor版——`A ⊗[ℚ] ℝ` 上の
`StandardEtalePair` の `cond` を `A ⊗[ℚ] R.1` へ降ろす。単射性の根拠が
(元の版の `Polynomial R` の自由性ではなく)`tensor_map_injective_of_flat`
(`ℚ`-平坦性)になる点だけが違う。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_fg_subalgebra_tensor_standardEtaleCond (A : Type) [CommRing A] [Algebra ℚ A]
    (f g p₁ p₂ : Polynomial (A ⊗[ℚ] ℝ)) (n : ℕ)
    (hcond : Polynomial.derivative f * p₁ + f * p₂ = g ^ n) :
    ∃ (R : FgSubalgebra ℚ ℝ) (f₀ g₀ p₁₀ p₂₀ : Polynomial (A ⊗[ℚ] R.1)),
      f₀.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom = f ∧
      g₀.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom = g ∧
      p₁₀.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom = p₁ ∧
      p₂₀.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom = p₂ ∧
      Polynomial.derivative f₀ * p₁₀ + f₀ * p₂₀ = g₀ ^ n := by
  classical
  obtain ⟨R, q, hq⟩ := exists_fg_subalgebra_tensor_polynomial_family A
    (![f, g, p₁, p₂] : Fin 4 → Polynomial (A ⊗[ℚ] ℝ))
  have hf : (q 0).map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom
      = f := by simpa using hq 0
  have hg : (q 1).map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom
      = g := by simpa using hq 1
  have hp1 : (q 2).map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom
      = p₁ := by simpa using hq 2
  have hp2 : (q 3).map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom
      = p₂ := by simpa using hq 3
  refine ⟨R, q 0, q 1, q 2, q 3, hf, hg, hp1, hp2, ?_⟩
  have hinj : Function.Injective
      (Polynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom) :=
    Polynomial.map_injective _ (tensor_map_injective_of_flat A R)
  apply hinj
  rw [Polynomial.map_add, Polynomial.map_mul, Polynomial.map_mul, Polynomial.map_pow,
    ← Polynomial.derivative_map, hf, hg, hp1, hp2, hcond]

open scoped TensorProduct in
/-- **`StandardEtalePair` は `A ⊗[ℚ] ℝ` から `A ⊗[ℚ] R.1`(ある
`R : FgSubalgebra ℚ ℝ`)へ降ろせる**——`exists_fg_subalgebra_
standardEtalePair` のtensor版。`Lemma 4.1` の「1アフィン片の降下」構成の
核となる部品(`ExtLimit.lean` の `piecePullbackIso` と組み合わせて使う)。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_fg_subalgebra_tensor_standardEtalePair (A : Type) [CommRing A] [Algebra ℚ A]
    (P : StandardEtalePair (A ⊗[ℚ] ℝ)) :
    ∃ (R : FgSubalgebra ℚ ℝ) (P₀ : StandardEtalePair (A ⊗[ℚ] R.1)),
      P₀.f.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom = P.f ∧
      P₀.g.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom = P.g := by
  obtain ⟨p₁, p₂, n, hcond⟩ := P.cond
  obtain ⟨R, f₀, g₀, p₁₀, p₂₀, hf, hg, hp1, hp2, hcond₀⟩ :=
    exists_fg_subalgebra_tensor_standardEtaleCond A P.f P.g p₁ p₂ n hcond
  have hinj : Function.Injective
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom :=
    tensor_map_injective_of_flat A R
  have hmonic : f₀.Monic := by
    refine (hinj.monic_map_iff (p := f₀)).mpr ?_
    rw [hf]; exact P.monic_f
  exact ⟨R, ⟨f₀, hmonic, g₀, ⟨p₁₀, p₂₀, n, hcond₀⟩⟩, hf, hg⟩

open scoped TensorProduct in
/-- **`exists_fg_subalgebra_tensor_standardEtalePair` と
`standardEtalePairRingBaseChange` を組み合わせた完成形**——
`StandardEtalePair (A ⊗[ℚ] ℝ)` は、ある有限段階 `R` 上の
`StandardEtalePair (A ⊗[ℚ] R.1)` の**base change として文字通り
一致する**(`P₀.Ring` を `A ⊗[ℚ] ℝ` へ係数拡大すると `P.Ring` に同型)。
`Lemma 4.1` の「1アフィン片の降下」で、`c.α` の局所片(標準エタール表示)
を有限段階へ降ろし、かつ「降ろした先の base change が元の局所片に一致
する」ことまで保証する、核心の完成形。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_fg_subalgebra_tensor_standardEtalePair_baseChange (A : Type) [CommRing A]
    [Algebra ℚ A] (P : StandardEtalePair (A ⊗[ℚ] ℝ)) :
    ∃ (R : FgSubalgebra ℚ ℝ) (P₀ : StandardEtalePair (A ⊗[ℚ] R.1)),
      letI : Algebra (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ) :=
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom.toAlgebra
      Nonempty (TensorProduct (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ) P₀.Ring ≃ₐ[A ⊗[ℚ] ℝ] P.Ring) := by
  obtain ⟨R, P₀, hf, hg⟩ := exists_fg_subalgebra_tensor_standardEtalePair A P
  letI : Algebra (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ) :=
    (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom.toAlgebra
  refine ⟨R, P₀, ?_⟩
  have hPeq : P₀.map (algebraMap (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ)) = P := by
    have hf' : (P₀.map (algebraMap (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ))).f = P.f := by
      rw [StandardEtalePair.map_f]; exact hf
    have hg' : (P₀.map (algebraMap (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ))).g = P.g := by
      rw [StandardEtalePair.map_g]; exact hg
    cases hpm : P₀.map (algebraMap (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ)) with
    | mk f' monic_f' g' cond' =>
      cases P with
      | mk f monic_f g cond =>
        rw [hpm] at hf' hg'
        simp only at hf' hg'
        subst hf' hg'
        rfl
  exact ⟨hPeq ▸ standardEtalePairRingBaseChange (K := A ⊗[ℚ] ℝ) P₀⟩

/-- **`StandardEtalePair.map` は合成と可換**——`(P.map φ).map ψ = P.map
(ψ.comp φ)`。`f`・`g` の一致(`Polynomial.map_map`)から構造体そのものの
一致を出す(`exists_fg_subalgebra_tensor_standardEtalePair_baseChange`と
同じ技法)。`Lemma 4.1` の「複数の標準エタール片を横断した細分段階の
合流」で、ある段階 `R_i` で得た局所モデルを、より粗い共通段階 `R'`
(`R_i ≤ R'`)へ**移送しても base change が変わらない**ことを示すのに使う。

★**sorry 無し**。標準3公理のみ。 -/
theorem StandardEtalePair.map_map {R S T : Type} [CommRing R] [CommRing S] [CommRing T]
    (φ : R →+* S) (ψ : S →+* T) (P : StandardEtalePair R) :
    (P.map φ).map ψ = P.map (ψ.comp φ) := by
  have hf : ((P.map φ).map ψ).f = (P.map (ψ.comp φ)).f := by
    rw [StandardEtalePair.map_f, StandardEtalePair.map_f, StandardEtalePair.map_f,
      Polynomial.map_map]
  have hg : ((P.map φ).map ψ).g = (P.map (ψ.comp φ)).g := by
    rw [StandardEtalePair.map_g, StandardEtalePair.map_g, StandardEtalePair.map_g,
      Polynomial.map_map]
  cases hL : (P.map φ).map ψ with
  | mk f' monic_f' g' cond' =>
    cases hR : P.map (ψ.comp φ) with
    | mk f monic_f g cond =>
      rw [hL] at hf hg
      rw [hR] at hf hg
      simp only at hf hg
      subst hf hg
      rfl

open scoped TensorProduct in
/-- **有限段階 `R` 上の降下 `P₀` は、より粗い共通段階 `R'`(`R ≤ R'`)へ
移送しても base change が変わらない**——`StandardEtalePair.map_map`
(合成との可換性)+`Algebra.TensorProduct.map_comp`(mathlib、`⊗`の関手性)
で `algebraMap(A⊗R.1)(A⊗ℝ) = algebraMap(A⊗R'.1)(A⊗ℝ) ∘ (A⊗inclusion)`
を示すだけ。`Lemma 4.1`の「複数の標準エタール片を横断した細分段階の合流」
——`exists_fgSubalgebra_upperBound`で得た共通段階`R'`へ、個々の片の
降下先`R_i`をこの補題で移送し、すべて同じ`R'`の上で扱えるようにする。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_fg_subalgebra_tensor_standardEtalePair_promote (A : Type) [CommRing A]
    [Algebra ℚ A] (R R' : FgSubalgebra ℚ ℝ) (hle : R ≤ R') (P₀ : StandardEtalePair (A ⊗[ℚ] R.1))
    (P : StandardEtalePair (A ⊗[ℚ] ℝ))
    (hP : P₀.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom = P) :
    (P₀.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hle)).toRingHom).map
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom = P := by
  rw [StandardEtalePair.map_map]
  rw [← hP]
  have hval : (Subalgebra.val R'.1).comp (Subalgebra.inclusion hle) = Subalgebra.val R.1 := rfl
  have hid : (AlgHom.id ℚ A).comp (AlgHom.id ℚ A) = AlgHom.id ℚ A := by ext; simp
  have hcomp := Algebra.TensorProduct.map_comp (AlgHom.id ℚ A) (AlgHom.id ℚ A)
    (Subalgebra.val R'.1) (Subalgebra.inclusion hle)
  rw [hid, hval] at hcomp
  have hcomp' : (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom.comp
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hle)).toRingHom =
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom :=
    congrArg AlgHom.toRingHom hcomp.symm
  rw [← hcomp']

/-! ## 重なりの比較(作業単位1)へ向けた第一歩——`Away(x*y)` の特徴づけ

`Lemma 4.1`の「比較射の構成」(`corrhyp-goal.md` §4のロードマップ、
作業単位1)最初の代数的補題:`D(x)∩D(y) = D(xy)`を、
`Localization.Away(x*y)`が(`Submonoid.powers(xy)`だけでなく)
`Submonoid.closure{x,y}`に関しても局所化になっている、という形で
特徴づける——`Away(x)`をさらに`y`で局所化したものとの比較の土台。 -/

/-- `Submonoid.powers(x*y) ≤ Submonoid.closure{x,y}`——`(xy)^n = x^n·y^n`
は明らかに`x,y`から生成される部分モノイドに入る。 -/
theorem powers_le_closure {R : Type} [CommRing R] (x y : R) :
    Submonoid.powers (x * y) ≤ (Submonoid.closure ({x, y} : Set R)) := by
  rw [Submonoid.powers_le]
  exact Submonoid.mul_mem _ (Submonoid.subset_closure (by simp)) (Submonoid.subset_closure (by simp))

/-- `closure{x,y}`の任意の元`z`について、ある`m`で`m*z`が`powers(xy)`に
入る——`closure_induction`で`x`・`y`・`1`・積の4ケースに分解するだけ。 -/
theorem exists_mul_mem_powers_of_closure {R : Type} [CommRing R] (x y : R) (z : R)
    (hz : z ∈ Submonoid.closure ({x, y} : Set R)) :
    ∃ m : R, m * z ∈ Submonoid.powers (x * y) := by
  induction hz using Submonoid.closure_induction with
  | mem a ha =>
      rcases ha with rfl | rfl
      · exact ⟨y, 1, by ring⟩
      · exact ⟨x, 1, by ring⟩
  | one => exact ⟨1, 0, by ring⟩
  | mul a b ha hb iha ihb =>
      obtain ⟨ma, na, hna⟩ := iha
      obtain ⟨mb, nb, hnb⟩ := ihb
      simp only at hna hnb
      refine ⟨ma * mb, na + nb, ?_⟩
      show (x*y) ^ (na + nb) = ma * mb * (a * b)
      rw [pow_add, hna, hnb]
      ring

/-- **`Localization.Away(x*y)` は `Submonoid.closure{x,y}` に関しても
局所化になっている**——`powers(xy)`と`closure{x,y}`が同じ飽和
(saturation)を持つことから(`IsLocalization.isLocalization_of_is_
exists_mul_mem`)。`Away(x)`をさらに`y`で局所化したものとの比較の
第一歩(作業単位1)。

★**sorry 無し**。標準3公理のみ。 -/
theorem isLocalization_closure_away_mul {R : Type} [CommRing R] (x y : R) :
    IsLocalization (Submonoid.closure ({x, y} : Set R)) (Localization.Away (x * y)) :=
  IsLocalization.isLocalization_of_is_exists_mul_mem (Localization.Away (x*y))
    (Submonoid.powers (x*y)) (Submonoid.closure ({x, y} : Set R))
    (powers_le_closure x y)
    (fun z => exists_mul_mem_powers_of_closure x y z.1 z.2)

/-! ## `StandardEtalePair.Ring` の元を有限段階へ降ろす——作業単位1(b)の土台

作業単位1(比較射の構成)で本当に要るのは「捩れの場合分け」ではなく、
「重なりを指定する元(比較射の分母)自体を有限段階`R'`の元として認識
すること」だと判明した(`corrhyp-goal.md`参照)。その第一歩として、
`P.Ring`(`P : StandardEtalePair (A⊗[ℚ]ℝ)`、`Bivariate` 多項式環
`(A⊗ℝ)[X][Y]` の商)の元を有限個同時に有限段階へ降ろす道具を用意する
——`exists_fg_subalgebra_tensor_finset`の2変数多項式版。 -/

open scoped TensorProduct in
/-- `exists_fg_subalgebra_tensor_finset` の2変数多項式版——`Polynomial
(Polynomial (A⊗[ℚ]ℝ))`(`StandardEtalePair.Ring`が経由する`Bivariate`
多項式環そのもの)の有限個の元は、ある共通の`R : FgSubalgebra ℚ ℝ`から
作った`Polynomial (Polynomial (A⊗[ℚ]R.1))`の元の像として同時に書ける。
係数(各多項式の係数のさらに係数)を1つの`Finset (A⊗[ℚ]ℝ)`へ平坦化して
`exists_fg_subalgebra_tensor_finset`に渡し、得られた降下先の係数から
2重の`monomial`和で元の多項式を組み直す——`exists_fg_subalgebra_tensor_
polynomial_family`と同じ手筋を1段深くしただけ。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_fg_subalgebra_tensor_bivariate_finset (A : Type) [CommRing A] [Algebra ℚ A]
    (s : Finset (Polynomial (Polynomial (A ⊗[ℚ] ℝ)))) :
    ∃ R : FgSubalgebra ℚ ℝ, ∀ q ∈ s, ∃ q₀ : Polynomial (Polynomial (A ⊗[ℚ] R.1)),
      q₀.map (Polynomial.mapRingHom
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom) = q := by
  classical
  set allLeaves : Finset (A ⊗[ℚ] ℝ) :=
    s.biUnion (fun q => q.support.biUnion
      (fun i => (q.coeff i).support.image (fun j => (q.coeff i).coeff j))) with hallLeaves
  obtain ⟨R, hR⟩ := exists_fg_subalgebra_tensor_finset A allLeaves
  refine ⟨R, ?_⟩
  intro q hq
  set φ := (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom with hφ
  choose c hc using fun t : Σ i : q.support, (q.coeff (i:ℕ)).support =>
    hR ((q.coeff (t.1:ℕ)).coeff (t.2:ℕ)) (by
      rw [hallLeaves]
      refine Finset.mem_biUnion.mpr ⟨q, hq, ?_⟩
      refine Finset.mem_biUnion.mpr ⟨(t.1 : ℕ), t.1.2, ?_⟩
      exact Finset.mem_image_of_mem _ t.2.2)
  refine ⟨q.support.attach.sum (fun i => Polynomial.monomial (i : ℕ)
    ((q.coeff (i:ℕ)).support.attach.sum (fun j =>
      Polynomial.monomial (j : ℕ) (c ⟨i, j⟩)))), ?_⟩
  have key : ∀ i : q.support, (Polynomial.monomial (i : ℕ)
      ((q.coeff (i:ℕ)).support.attach.sum (fun j => Polynomial.monomial (j : ℕ) (c ⟨i, j⟩)))).map
      (Polynomial.mapRingHom φ) =
      Polynomial.monomial (i : ℕ) (q.coeff (i : ℕ)) := by
    intro i
    rw [Polynomial.map_monomial]
    have hcoe : (⇑(Polynomial.mapRingHom φ) : Polynomial (A ⊗[ℚ] R.1) → Polynomial (A ⊗[ℚ] ℝ))
        = Polynomial.map φ := rfl
    rw [hcoe]
    congr 1
    rw [Polynomial.map_sum]
    have key2 : ∀ j : (q.coeff (i:ℕ)).support,
        (Polynomial.monomial (j : ℕ) (c ⟨i, j⟩)).map φ =
        Polynomial.monomial (j : ℕ) ((q.coeff (i:ℕ)).coeff (j : ℕ)) := by
      intro j
      rw [Polynomial.map_monomial]
      exact congrArg (Polynomial.monomial (j : ℕ)) (hc ⟨i, j⟩)
    rw [Finset.sum_congr rfl (fun j (_ : j ∈ (q.coeff (i:ℕ)).support.attach) => key2 j)]
    rw [Finset.sum_attach (q.coeff (i:ℕ)).support
      (fun j => Polynomial.monomial j ((q.coeff (i:ℕ)).coeff j))]
    exact Polynomial.sum_monomial_eq (q.coeff (i:ℕ))
  rw [Polynomial.map_sum]
  rw [Finset.sum_congr rfl (fun i (_ : i ∈ q.support.attach) => key i)]
  rw [Finset.sum_attach q.support (fun i => Polynomial.monomial i (q.coeff i))]
  exact Polynomial.sum_monomial_eq q

open scoped TensorProduct in
/-- **`exists_fg_subalgebra_tensor_bivariate_finset`の`MvPolynomial ι`
版**——`ι`変数多項式(有限個の族)を、その係数(有限個、`support`+`coeff`
で取り出せる)を`exists_fg_subalgebra_tensor_finset`で共通の`R`へ
降ろしてから`monomial`の和で再構成することで、単一の`R`へ降ろす。
`Algebra.FinitePresentation`が与える`MvPolynomial`商としての表示
(`Algebra.FinitePresentation.iff_quotient_mvPolynomial'`、mathlib)を
`R`レベルへ降ろす計画(`Γ(C,piece)`のような有限表示な代数を、個別の
元ごとではなく一度に降ろす、2026-09-04の設計再考)の核心部品。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_fg_subalgebra_tensor_mvPolynomial_finset (A : Type) [CommRing A] [Algebra ℚ A]
    {ι : Type} (s : Finset (MvPolynomial ι (A ⊗[ℚ] ℝ))) :
    ∃ R : FgSubalgebra ℚ ℝ, ∀ q ∈ s, ∃ q₀ : MvPolynomial ι (A ⊗[ℚ] R.1),
      MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom q₀ = q := by
  classical
  set allCoeffs : Finset (A ⊗[ℚ] ℝ) :=
    s.biUnion (fun q => q.support.image (fun m => q.coeff m)) with hallCoeffs
  obtain ⟨R, hR⟩ := exists_fg_subalgebra_tensor_finset A allCoeffs
  refine ⟨R, ?_⟩
  intro q hq
  choose c hc using fun m : q.support => hR (q.coeff (m : ι →₀ ℕ)) (by
    rw [hallCoeffs]
    exact Finset.mem_biUnion.mpr ⟨q, hq, Finset.mem_image_of_mem _ m.2⟩)
  refine ⟨q.support.attach.sum (fun m => MvPolynomial.monomial (m : ι →₀ ℕ) (c m)), ?_⟩
  rw [map_sum]
  have key : ∀ m : q.support,
      MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom
        (MvPolynomial.monomial (m : ι →₀ ℕ) (c m)) =
      MvPolynomial.monomial (m : ι →₀ ℕ) (q.coeff (m : ι →₀ ℕ)) := by
    intro m
    rw [MvPolynomial.map_monomial]
    show MvPolynomial.monomial (m : ι →₀ ℕ) ((Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)) (c m))
      = _
    rw [hc m]
  rw [Finset.sum_congr rfl (fun m _ => key m), Finset.sum_attach q.support (fun m => MvPolynomial.monomial m (q.coeff m))]
  exact (q.support_sum_monomial_coeff)

open scoped TensorProduct in
/-- `MvPolynomial ι R ⊗[R] A ≃+* MvPolynomial ι A`の橋渡し(補助部品)
——`Algebra.TensorProduct.comm`(係数順序を`A ⊗[R] MvPolynomial ι R`へ
入れ替える)と`MvPolynomial.algebraTensorAlgEquiv`(mathlib、`A⊗
MvPolynomial σ R ≃ₐ[A] MvPolynomial σ A`)の合成。`quotient_mvPolynomial_
baseChange`の中核部品。 -/
noncomputable def quotientMvPolynomialBaseChangeRingEquivAux (R A : Type) [CommRing R] [CommRing A] [Algebra R A]
    (ι : Type) : MvPolynomial ι R ⊗[R] A ≃+* MvPolynomial ι A :=
  (Algebra.TensorProduct.comm R (MvPolynomial ι R) A).toRingEquiv.trans
    (MvPolynomial.algebraTensorAlgEquiv R A).toRingEquiv

open scoped TensorProduct in
/-- `quotientMvPolynomialBaseChangeRingEquivAux`が、`MvPolynomial ι R`
の`algebraMap`による埋め込み(`x ↦ x⊗ₜ1`)を`MvPolynomial.map (algebraMap
R A)`と両立させること——`Algebra.TensorProduct.algebraMap_apply`+
`Algebra.TensorProduct.comm_tmul`+`MvPolynomial.algebraTensorAlgEquiv_
symm_map`(いずれもmathlib)を繋ぐだけ。 -/
theorem quotientMvPolynomialBaseChangeRingEquivAux_comp_algebraMap (R A : Type) [CommRing R] [CommRing A]
    [Algebra R A] (ι : Type) (x : MvPolynomial ι R) :
    quotientMvPolynomialBaseChangeRingEquivAux R A ι
      (algebraMap (MvPolynomial ι R) (MvPolynomial ι R ⊗[R] A) x) = MvPolynomial.map (algebraMap R A) x := by
  unfold quotientMvPolynomialBaseChangeRingEquivAux
  rw [Algebra.TensorProduct.algebraMap_apply]
  show (MvPolynomial.algebraTensorAlgEquiv R A)
    ((Algebra.TensorProduct.comm R (MvPolynomial ι R) A) (x ⊗ₜ[R] 1)) = _
  rw [Algebra.TensorProduct.comm_tmul]
  exact (MvPolynomial.algebraTensorAlgEquiv_symm_map R A x) ▸ (AlgEquiv.apply_symm_apply _ _)

open scoped TensorProduct in
/-- **`MvPolynomial`の商とbase changeが可換であること**——`(MvPolynomial
ι R ⧸ I) ⊗[R] A ≃+* MvPolynomial ι A ⧸ Ideal.map (MvPolynomial.map
(algebraMap R A)) I`。`Algebra.TensorProduct.quotientTensorEquiv`
(mathlib、商とテンソル積の可換性)を`quotientMvPolynomialBaseChange
RingEquivAux`で`MvPolynomial ι A`の形へ書き換えたもの。`Γ(C,piece)`の
`R`レベルモデル`S_0`(`pieceAlgebra_R_model`、`ExtLimit.lean`)のbase
changeが実際に`Γ(C,piece)`を復元することを示す核心部品。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def quotient_mvPolynomial_baseChange (R A : Type) [CommRing R] [CommRing A] [Algebra R A] (ι : Type)
    (I : Ideal (MvPolynomial ι R)) :
    (MvPolynomial ι R ⧸ I) ⊗[R] A ≃+* MvPolynomial ι A ⧸ Ideal.map (MvPolynomial.map (algebraMap R A)) I := by
  refine (Algebra.TensorProduct.quotientTensorEquiv (R := R) R A (MvPolynomial ι R) I).toRingEquiv.trans ?_
  refine Ideal.quotientEquiv _ _ (quotientMvPolynomialBaseChangeRingEquivAux R A ι) ?_
  rw [Ideal.map_map]
  exact congrArg (Ideal.map · I)
    (RingHom.ext (fun x => (quotientMvPolynomialBaseChangeRingEquivAux_comp_algebraMap R A ι x).symm))

/-! ## `A ⊗[ℚ] R'.1 → A ⊗[ℚ] ℝ` は単射——「項目(d)の第二段」
(遷移データのRレベル降下)の鍵

`ℚ` が体なので任意の`ℚ`-加群(とくに`A`)は`Module.Flat`——`R'.1 ↪ ℝ`
という単射を`A`で(左)テンソルしても単射性が保たれる
(`Module.Flat.lTensor_preserves_injective_linearMap`、mathlib)。これは
`corrhyp-goal.md`2026-09-05で見積もった`RingHom.EssFiniteType`の
filtered colimit descent(`exists_comp_map_eq_of_isColimit`等)よりも
**直接的**な道筋——「ℝレベルで等しいなら、共通の精密化を探すまでもなく
**その場の`R'`で既に等しい**」という強い主張が、体上のベクトル空間の
平坦性から即座に出る。遷移同型の**一意性**(貼り合わせの整合性)は、
このまま`MvPolynomial`の係数ごとの単射性に直結する——2つの候補写像が
ℝへbase changeすると一致するなら、生成元の像(有限個)を係数として
持つ多項式として、その場で一致することが言える。 -/

open scoped TensorProduct in
/-- **`A ⊗[ℚ] R'.1 → A ⊗[ℚ] ℝ` は単射**——`ℚ`が体であること
(`Module.Flat ℚ A`が自動)+`Subalgebra.val R'.1`の単射性+
`Module.Flat.lTensor_preserves_injective_linearMap`(mathlib)。

★**sorry 無し**。標準3公理のみ。 -/
theorem algebraTensorMap_val_injective (A : Type) [CommRing A] [Algebra ℚ A] (R' : FgSubalgebra ℚ ℝ) :
    Function.Injective (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)) := by
  have hinj : Function.Injective (Subalgebra.val R'.1) := Subtype.val_injective
  intro x y h
  apply Module.Flat.lTensor_preserves_injective_linearMap (R := ℚ) (M := A)
    (Subalgebra.val R'.1).toLinearMap hinj
  show (LinearMap.lTensor A (Subalgebra.val R'.1).toLinearMap) x =
    (LinearMap.lTensor A (Subalgebra.val R'.1).toLinearMap) y
  rw [show (LinearMap.lTensor A (Subalgebra.val R'.1).toLinearMap) =
    (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toLinearMap from rfl]
  exact h

open scoped TensorProduct in
/-- `algebraTensorMap_val_injective`の`MvPolynomial`版——
`MvPolynomial.map`(mathlib)は単射な係数写像を単射に保つ。

★**sorry 無し**。標準3公理のみ。 -/
theorem mvPolynomial_map_algebraTensorMap_val_injective (A : Type) [CommRing A] [Algebra ℚ A]
    (R' : FgSubalgebra ℚ ℝ) (ι : Type) :
    Function.Injective (MvPolynomial.map (σ := ι)
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom) :=
  MvPolynomial.map_injective _ (algebraTensorMap_val_injective A R')

open scoped TensorProduct in
/-- **`R'`レベルの多項式が`ℝ`へbase changeして`0`になれば、`R'`レベルで
既に`0`だった**——`mvPolynomial_map_algebraTensorMap_val_injective`の
言い換え。遷移同型の候補(`R'`レベル)が、生成元の関係式をℝレベルで
満たすことが分かれば、それだけで`R'`レベルでも関係式を満たす
(well-definedness)ことを保証する、貼り合わせの核心部品。 -/
theorem mvPolynomial_algebraTensorMap_val_eq_zero_of_map_eq_zero (A : Type) [CommRing A] [Algebra ℚ A]
    (R' : FgSubalgebra ℚ ℝ) (ι : Type) (p : MvPolynomial ι (A ⊗[ℚ] R'.1))
    (h : MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom p = 0) :
    p = 0 :=
  (mvPolynomial_map_algebraTensorMap_val_injective A R' ι).eq_iff.mp (h.trans (map_zero _).symm)

open scoped TensorProduct in
/-- **`MvPolynomial.map`と`MvPolynomial.aeval`が可換であること**——
`ψ := algebraMap (A⊗R'.1)(A⊗ℝ)`による係数拡大と、生成元への代入
(`aeval`)の順序を入れ替えても同じ結果になる。`MvPolynomial.induction_on`
での3ケース(定数・和・単項式の積)による直接計算——mathlibに完成品の
1つの補題として無かったので手で組んだ。「候補の遷移写像(`R'`レベルの
`ev`)が、ℝへbase changeした後に関係式を満たす」ことと「`R'`レベルで
関係式を満たす」ことを結ぶ橋渡し。

★**sorry 無し**。標準3公理のみ。 -/
theorem mvPolynomial_map_aeval_comm (A : Type) [CommRing A] [Algebra ℚ A] (R' : FgSubalgebra ℚ ℝ) (ι ι' : Type)
    (p : MvPolynomial ι (A ⊗[ℚ] R'.1)) (ev : ι → MvPolynomial ι' (A ⊗[ℚ] R'.1)) :
    MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom
        (MvPolynomial.aeval ev p) =
      MvPolynomial.aeval (fun i => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom (ev i))
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom p) := by
  induction p using MvPolynomial.induction_on with
  | C a =>
    simp only [MvPolynomial.aeval_C, MvPolynomial.algebraMap_eq, MvPolynomial.map_C]
  | add p q hp hq =>
    rw [map_add, map_add, map_add, hp, hq, map_add]
  | mul_X p i hp =>
    simp only [map_mul, MvPolynomial.aeval_X, MvPolynomial.map_X]
    rw [hp]

open scoped TensorProduct in
/-- **遷移写像の候補(`R'`レベルの`ev`)がℝレベルで関係式`p`を満たすなら、
`R'`レベルで既に満たしている**——`mvPolynomial_map_aeval_comm`(係数拡大
とaevalの可換性)+`mvPolynomial_algebraTensorMap_val_eq_zero_of_map_eq_zero`
(単射性)を合成するだけ。「項目(d)の第二段」(遷移データのRレベル降下)
で実際に使う核心の補題——`ev`(候補の生成元の像、`exists_fg_subalgebra_
tensor_mvPolynomial_finset`で降ろしたもの)が、既知のℝレベルの遷移
写像が満たす関係式を、ℝへbase changeするまでもなく`R'`レベルで
すでに満たすことを保証する(well-definedness)。

★**sorry 無し**。標準3公理のみ。 -/
theorem mvPolynomial_aeval_eq_zero_of_map_aeval_eq_zero (A : Type) [CommRing A] [Algebra ℚ A]
    (R' : FgSubalgebra ℚ ℝ) (ι ι' : Type) (p : MvPolynomial ι (A ⊗[ℚ] R'.1))
    (ev : ι → MvPolynomial ι' (A ⊗[ℚ] R'.1))
    (h : MvPolynomial.aeval (fun i => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom (ev i))
      (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom p) = 0) :
    MvPolynomial.aeval ev p = 0 := by
  apply mvPolynomial_algebraTensorMap_val_eq_zero_of_map_eq_zero
  rw [mvPolynomial_map_aeval_comm]
  exact h

open scoped TensorProduct in
/-- **`R₀→R'`への昇格と、両方を`ℝ`へ送る写像が両立すること**——
`(A⊗R'.1へ)∘(R₀→R'への昇格) = (A⊗R₀.1からℝへ)`。`Algebra.TensorProduct.
map_comp`(mathlib、テンソル関手性)+`Subalgebra.val_comp_inclusion`
(mathlib、`T.val∘(inclusion h)=S.val`)を繋ぐだけ。`R'`レベルへ昇格した
元が、実際に元の`R₀`レベルでの値をℝへ送ったものと一致することを保証
する、貼り合わせの配線でくり返し使う基本事実。

★**sorry 無し**。標準3公理のみ。 -/
theorem algebraTensorMap_val_comp_inclusion (A : Type) [CommRing A] [Algebra ℚ A] {R₀ R' : FgSubalgebra ℚ ℝ}
    (h : R₀ ≤ R') :
    (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom.comp
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion h)).toRingHom =
    (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R₀.1)).toRingHom := by
  apply RingHom.ext
  intro x
  show (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1))
      ((Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion h)) x) = _
  rw [← AlgHom.comp_apply, ← Algebra.TensorProduct.map_comp, Subalgebra.val_comp_inclusion, AlgHom.id_comp]
  rfl

open scoped TensorProduct in
/-- **「項目(d)の第二段」の核心——ℝレベルで定義された生成元の対応
(有限個)が関係式を満たすなら、それは`R`レベルの候補写像の`R'`レベル
への昇格として実際に実現される**。`q:κ→MvPolynomial ι(A⊗R.1)`(片1の
関係式、`R`レベル)と`ψ:ι→MvPolynomial ι'(A⊗ℝ)`(既知のℝレベルの
遷移写像が生成元に送る値)から出発し、
1. `exists_fg_subalgebra_tensor_mvPolynomial_finset`で`ψ`の**有限個**の
   値(`ι`が`Fintype`)を共通の`R₀`へ降ろす(`ev₀`)。
2. `exists_fgSubalgebra_upperBound2`で`R`と`R₀`を共通の`R'`へ合流
   させる。
3. `ev₀`を`R'`へ昇格した`ev`が、実際に`ψ`を再現すること
   (`algebraTensorMap_val_comp_inclusion`)と、`q`(`R'`へ昇格したもの)
   の関係式を`R'`レベルで満たすこと(`mvPolynomial_aeval_eq_zero_of_
   map_aeval_eq_zero`、単射性経由)の両方を確認する。

これで、「ℝレベルで分かっている遷移写像の生成元の対応」から「実際に
`R'`レベルの候補写像」を**構成的に**取り出せる——`transitionElem`/
`gdT`/`cocycle`のRレベル版に相当する遷移データの降下が、個別の
GlueDataエンジニアリングの再構築ではなく、この1つの補題への
`Presentation`データの specialize として実現できる見込みが立った。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_mvPolynomial_quotient_ringHom_descend (A : Type) [CommRing A] [Algebra ℚ A]
    (R : FgSubalgebra ℚ ℝ) {ι κ ι' : Type} [Fintype ι]
    (q : κ → MvPolynomial ι (A ⊗[ℚ] R.1)) (ψ : ι → MvPolynomial ι' (A ⊗[ℚ] ℝ))
    (hψ : ∀ k, MvPolynomial.aeval ψ
      (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom (q k)) = 0) :
    ∃ (R' : FgSubalgebra ℚ ℝ) (hR : R ≤ R') (ev : ι → MvPolynomial ι' (A ⊗[ℚ] R'.1)),
      (∀ i, MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom (ev i)
        = ψ i) ∧
      (∀ k, MvPolynomial.aeval ev (MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR)).toRingHom (q k)) = 0) := by
  classical
  obtain ⟨R₀, hR₀spec⟩ := exists_fg_subalgebra_tensor_mvPolynomial_finset A (Finset.image ψ Finset.univ)
  choose ev₀ hev₀ using fun i : ι => hR₀spec (ψ i) (Finset.mem_image_of_mem _ (Finset.mem_univ i))
  obtain ⟨R', hR, hR₀⟩ := exists_fgSubalgebra_upperBound2 R R₀
  have hcomm : ∀ i, MvPolynomial.map
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom
      (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR₀)).toRingHom
        (ev₀ i)) = ψ i := by
    intro i
    rw [MvPolynomial.map_map, algebraTensorMap_val_comp_inclusion]
    exact hev₀ i
  refine ⟨R', hR, fun i => MvPolynomial.map
    (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR₀)).toRingHom (ev₀ i), hcomm, ?_⟩
  intro k
  apply mvPolynomial_aeval_eq_zero_of_map_aeval_eq_zero
  rw [show (fun i => MvPolynomial.map
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom
      (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR₀)).toRingHom
        (ev₀ i))) = ψ from funext hcomm]
  rw [MvPolynomial.map_map, algebraTensorMap_val_comp_inclusion]
  exact hψ k

/-! ## `StandardEtalePair.Ring` の元を有限段階へ降ろす——作業単位1(b)の完成

`exists_fg_subalgebra_tensor_bivariate_finset`を`StandardEtalePair.Ring`
(`Bivariate`多項式環の商)の実際の元(比較射の分母`g_l`等)へ接続する。
`equivMvPolynomialQuotient`(`MvPolynomial`経由)を使う迂回は不要——
`P.Ring`は`.refl`で生の`Bivariate`商と同一視できるので、`Ideal.
quotientMap`を直接使って`P₀.Ring →+* (P₀.map f).Ring`という環準同型を
組み立てるほうが単純。 -/

open scoped TensorProduct in
/-- `exists_fg_subalgebra_tensor_standardEtalePair`の**構造体そのものの
一致版**——`exists_fg_subalgebra_tensor_standardEtalePair_baseChange`の
証明内部で使っていた技法(`.f`・`.g`の一致→`cases`+`subst`+`rfl`)を
単独の補題として取り出した。`_promote`はこの形の入力を要求する。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_fg_subalgebra_tensor_standardEtalePair_mapEq (A : Type) [CommRing A]
    [Algebra ℚ A] (P : StandardEtalePair (A ⊗[ℚ] ℝ)) :
    ∃ (R : FgSubalgebra ℚ ℝ) (P₀ : StandardEtalePair (A ⊗[ℚ] R.1)),
      P₀.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom = P := by
  obtain ⟨R, P₀, hf, hg⟩ := exists_fg_subalgebra_tensor_standardEtalePair A P
  refine ⟨R, P₀, ?_⟩
  have hf' : (P₀.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom).f
      = P.f := by rw [StandardEtalePair.map_f]; exact hf
  have hg' : (P₀.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom).g
      = P.g := by rw [StandardEtalePair.map_g]; exact hg
  cases hpm : P₀.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom with
  | mk f' monic_f' g' cond' =>
    cases P with
    | mk f monic_f g cond =>
      rw [hpm] at hf' hg'
      simp only at hf' hg'
      subst hf' hg'
      rfl

open Polynomial in
/-- **`P.Ring →+* (P.map f).Ring`という自然な環準同型**——`P.Ring`は
`equivPolynomialQuotient`が`.refl`であることが示す通り生の`Bivariate`
多項式環の商そのものなので、`Ideal.quotientMap`に`Polynomial.mapRingHom
(Polynomial.mapRingHom f)`(2変数への持ち上げ)を渡すだけで組み立てられる
——`(P.map f).f = P.f.map f`等(`StandardEtalePair.map`の定義)から、
この写像がイデアルの生成元をちょうど対応する生成元へ送ることを確認する。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def standardEtalePairMapRingHom {R S : Type} [CommRing R] [CommRing S] (f : R →+* S)
    (P : StandardEtalePair R) : P.Ring →+* (P.map f).Ring :=
  Ideal.quotientMap _ (Polynomial.mapRingHom (Polynomial.mapRingHom f)) (by
    rw [Ideal.span_le]
    intro x hx
    simp only [Set.mem_insert_iff, Set.mem_singleton_iff] at hx
    rw [SetLike.mem_coe, Ideal.mem_comap]
    rcases hx with rfl | rfl
    · show Polynomial.mapRingHom (Polynomial.mapRingHom f) (Polynomial.C P.f) ∈ _
      rw [show (Polynomial.mapRingHom (Polynomial.mapRingHom f)) (Polynomial.C P.f)
        = Polynomial.C (P.f.map f) from Polynomial.map_C _]
      exact Ideal.subset_span (by simp [StandardEtalePair.map_f])
    · show Polynomial.mapRingHom (Polynomial.mapRingHom f)
        (Polynomial.X * Polynomial.C P.g - 1) ∈ _
      rw [map_sub, map_mul, map_one]
      rw [show (Polynomial.mapRingHom (Polynomial.mapRingHom f)) Polynomial.X = Polynomial.X from
        Polynomial.map_X _]
      rw [show (Polynomial.mapRingHom (Polynomial.mapRingHom f)) (Polynomial.C P.g)
        = Polynomial.C (P.g.map f) from Polynomial.map_C _]
      exact Ideal.subset_span (by simp [StandardEtalePair.map_g]))

/-- `standardEtalePairMapRingHom`は`Ideal.Quotient.mk`と自然に可換
(`Ideal.quotientMap_mk`の直接の帰結、`FunLike`適用と`.map`dot記法の
不一致(lean-idioms #26)を`rfl`で吸収するだけ)。

★**sorry 無し**。標準3公理のみ。 -/
theorem standardEtalePairMapRingHom_mk {R S : Type} [CommRing R] [CommRing S] (f : R →+* S)
    (P : StandardEtalePair R) (q : Polynomial (Polynomial R)) :
    standardEtalePairMapRingHom f P (Ideal.Quotient.mk _ q) =
      Ideal.Quotient.mk _ (q.map (Polynomial.mapRingHom f)) := by
  have hcoe : (⇑(Polynomial.mapRingHom (Polynomial.mapRingHom f)) :
      Polynomial (Polynomial R) → Polynomial (Polynomial S)) =
      Polynomial.map (Polynomial.mapRingHom f) := rfl
  show Ideal.quotientMap _ _ _ (Ideal.Quotient.mk _ q) = _
  rw [Ideal.quotientMap_mk, hcoe]

open scoped TensorProduct in
/-- `Ideal (Polynomial (Polynomial (A ⊗[ℚ] ℝ)))`の`IsTwoSided`(可換環では
自動のはずの事実)が、`Polynomial`の2重ネスト+`A ⊗[ℚ] ℝ`という一般の
底環の組み合わせで通常の`inferInstance`では**候補の組み合わせ爆発**により
失敗する(`CommRing`インスタンスの探索が`Polynomial.commRing`を再帰的に
何度も試して失敗をキャッシュし、正しい経路へたどり着けない)ため、
`letI`で`CommRing`の連鎖を1段ずつ明示して与えてから`infer_instance`する
——この`instance`自体をグローバルに登録しておくと、以後この形の`Ideal`
を扱う任意の場所(`Ideal.Quotient.mk_surjective`等)で自動的に見つかる
ようになる。`tools/lean-idioms.md`に項目として追記予定。

★**sorry 無し**。標準3公理のみ。 -/
instance bivariateIsTwoSided {A : Type} [CommRing A] [Algebra ℚ A]
    (I : Ideal (Polynomial (Polynomial (A ⊗[ℚ] ℝ)))) : I.IsTwoSided := by
  letI i1 : CommRing (A ⊗[ℚ] ℝ) := inferInstance
  letI i2 : CommRing (Polynomial (A ⊗[ℚ] ℝ)) := @Polynomial.commRing (A ⊗[ℚ] ℝ) i1
  letI i3 : CommRing (Polynomial (Polynomial (A ⊗[ℚ] ℝ))) :=
    @Polynomial.commRing (Polynomial (A ⊗[ℚ] ℝ)) i2
  infer_instance

open scoped TensorProduct in
set_option maxHeartbeats 4000000 in
/-- **作業単位1(b)の完成——`P.Ring`の任意の元は有限段階へ降ろせる**。
`P : StandardEtalePair (A⊗[ℚ]ℝ)`の任意の元`z`について、ある有限段階
`R`・`P₀ : StandardEtalePair (A⊗[ℚ]R.1)`・`z₀ : P₀.Ring`が存在し、
`P₀`の base change が`P`に一致し(`_mapEq`)、かつ`z₀`の
`standardEtalePairMapRingHom`による像が`z`に一致する(`HEq`——`P₀.map f`
と`P`が値として等しいだけで型としては`▸`が必要になり、depth の深い
`Polynomial (Polynomial (A⊗ℚℝ))` 上で`▸`を文越しに使うと`whnf`が
combinatorial explosion で timeout するため、型の一致を要求しない`HEq`
で述べる)。手順: `z`を`Ideal.Quotient.mk`の代表元`q`へ持ち上げ
(`bivariateIsTwoSided`のおかげで`mk_surjective`が通る)、
`exists_fg_subalgebra_tensor_bivariate_finset`で`q`の係数を降ろし、
`P`自身の降下(`_mapEq`)と共通の上界`R'`(`⊔`)へ合流させ、
`standardEtalePairMapRingHom_mk`で像を計算する。

`Lemma 4.1`の「比較射の構成」(作業単位1)で、`f_l`・`f_m`(各局所片の
`StandardEtalePresentation`の生成元由来の元)自体を有限段階の元として
認識するのに直接使える。

★**sorry 無し**。標準3公理のみ。`lake build`で重い(約40秒、
`Ideal (Polynomial (Polynomial (A⊗ℚℝ)))`のインスタンス探索の深さのため
`maxHeartbeats`を上げている)が、有限時間で完全に閉じる。 -/
theorem exists_fg_subalgebra_tensor_standardEtale_elem (A : Type) [CommRing A] [Algebra ℚ A]
    (P : StandardEtalePair (A ⊗[ℚ] ℝ)) (z : P.Ring) :
    ∃ (R : FgSubalgebra ℚ ℝ) (P₀ : StandardEtalePair (A ⊗[ℚ] R.1)) (z₀ : P₀.Ring),
      P₀.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom = P ∧
      HEq (standardEtalePairMapRingHom
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom P₀ z₀) z := by
  obtain ⟨R₁, P₁, hP₁⟩ := exists_fg_subalgebra_tensor_standardEtalePair_mapEq A P
  obtain ⟨q, hq⟩ := Ideal.Quotient.mk_surjective z
  obtain ⟨R₂, hR₂⟩ := exists_fg_subalgebra_tensor_bivariate_finset A ({q} : Finset _)
  obtain ⟨q₀, hq₀⟩ := hR₂ q (Finset.mem_singleton_self q)
  set R' : FgSubalgebra ℚ ℝ := ⟨R₁.1 ⊔ R₂.1, fg_sup R₁.1 R₂.1 R₁.2 R₂.2⟩ with hR'def
  have hle1 : R₁ ≤ R' := show (R₁.1 : Subalgebra ℚ ℝ) ≤ R'.1 from le_sup_left
  have hle2 : R₂ ≤ R' := show (R₂.1 : Subalgebra ℚ ℝ) ≤ R'.1 from le_sup_right
  have hP₁' := exists_fg_subalgebra_tensor_standardEtalePair_promote A R₁ R' hle1 P₁ P hP₁
  refine ⟨R', P₁.map (Algebra.TensorProduct.map (AlgHom.id ℚ A)
    (Subalgebra.inclusion hle1)).toRingHom,
    Ideal.Quotient.mk _ (q₀.map (Polynomial.mapRingHom
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hle2)).toRingHom)),
    hP₁', ?_⟩
  rw [standardEtalePairMapRingHom_mk]
  have hval : (Subalgebra.val R'.1).comp (Subalgebra.inclusion hle2) = Subalgebra.val R₂.1 := rfl
  have hid : (AlgHom.id ℚ A).comp (AlgHom.id ℚ A) = AlgHom.id ℚ A := by ext; simp
  have hcomp := Algebra.TensorProduct.map_comp (AlgHom.id ℚ A) (AlgHom.id ℚ A)
    (Subalgebra.val R'.1) (Subalgebra.inclusion hle2)
  rw [hid, hval] at hcomp
  have hcomp' : (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom.comp
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hle2)).toRingHom =
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R₂.1)).toRingHom :=
    congrArg AlgHom.toRingHom hcomp.symm
  have hmm : (q₀.map (Polynomial.mapRingHom
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hle2)).toRingHom)).map
      (Polynomial.mapRingHom
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom) =
      q₀.map (Polynomial.mapRingHom
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R₂.1)).toRingHom) := by
    rw [Polynomial.map_map, Polynomial.mapRingHom_comp, hcomp']
  rw [hmm, hq₀]
  clear_value R'
  subst hP₁'
  rw [hq]

/-! ## `B`の元を`StandardEtalePresentation`のRingへ座標として運ぶ
——複数ピースの重なりの比較(work unit 3の続き)への橋渡し

複数の標準エタール片`D(f_i)`を貼り合わせるには、`f_j`(別の局所化元)が
`D(f_i)`(つまり`P_i.Ring`)の中でどんな元に対応するかを知る必要がある。
`StandardEtalePresentation.exists_mul_aeval_x_g_pow_eq_aeval_x`(mathlib、
「`S`の任意の元は`g(X)`のベキを掛ければ`X`の多項式になる」)を
`Pres.equivRing`で両辺に運ぶだけで、この「座標」を取り出せる。 -/

/-- **`StandardEtalePresentation`の元は、そのRing上で「`g`のベキ倍すれば
Xの多項式になる」という座標表示を持つ**——`exists_mul_aeval_x_g_pow_eq_
aeval_x`(mathlib)を`Pres.equivRing`(AlgHom)で両辺に運び、`aeval`の
自然性(`Polynomial.aeval_algHom_apply`)で`Pres.x`を`Pres.X`へ変換する
だけ。`D(f_i)`・`D(f_j)`の重なりを`P_i.Ring`の中の基本開集合として
記述するための第一歩——ここで得られる`p`(1変数多項式、`Polynomial R`)は
`exists_fg_subalgebra_tensor_polynomial_family`でそのまま有限段階へ
降ろせる。

★**sorry 無し**。標準3公理のみ。 -/
theorem standardEtalePresentation_exists_coord {R S : Type} [CommRing R] [CommRing S] [Algebra R S]
    (Pres : StandardEtalePresentation R S) (x : S) :
    ∃ (p : Polynomial R) (n : ℕ),
      Pres.equivRing x * (Polynomial.aeval Pres.P.X Pres.P.g) ^ n = Polynomial.aeval Pres.P.X p := by
  obtain ⟨p, n, hpn⟩ := Pres.exists_mul_aeval_x_g_pow_eq_aeval_x x
  refine ⟨p, n, ?_⟩
  have h2 := congrArg Pres.equivRing hpn
  rw [map_mul, map_pow, ← Polynomial.aeval_algHom_apply, ← Polynomial.aeval_algHom_apply,
    Pres.equivRing_x] at h2
  exact h2

/-- 単元の`PrimeSpectrum.basicOpen`は全体——mathlibに直接の名前が
見当たらなかったので自作(CorrHyp非依存の一般的事実)。 -/
theorem PrimeSpectrum.basicOpen_eq_top_of_isUnit {R : Type} [CommRing R] {u : R} (hu : IsUnit u) :
    PrimeSpectrum.basicOpen u = ⊤ := by
  apply TopologicalSpace.Opens.ext
  simp only [TopologicalSpace.Opens.coe_top]
  rw [Set.eq_univ_iff_forall]
  intro p
  rw [SetLike.mem_coe, PrimeSpectrum.mem_basicOpen]
  intro hmem
  exact p.2.ne_top (p.asIdeal.eq_top_of_isUnit_mem hmem hu)

/-- 単元倍しても`basicOpen`は変わらない——`basicOpen_mul`+
`basicOpen_eq_top_of_isUnit`から。重なりの遷移開集合を同定する
`standardEtalePresentation_transitionOpen_eq`の核心部品。 -/
theorem PrimeSpectrum.basicOpen_mul_isUnit {R : Type} [CommRing R] (a : R) {u : R} (hu : IsUnit u) :
    PrimeSpectrum.basicOpen (a * u) = PrimeSpectrum.basicOpen a := by
  rw [PrimeSpectrum.basicOpen_mul, PrimeSpectrum.basicOpen_eq_top_of_isUnit hu, inf_top_eq]

/-- **重なりの比較射に要る「遷移開集合」の同定**——`Pres`の元`x`
(たとえば`f_j`の`Localization.Away f_i`への像)を`Pres.equivRing`で
`Pres.P.Ring`へ運んだ像の基本開集合は、座標多項式`p`(`standardEtalePresentation_
exists_coord`で得た)を`X`で評価した元の基本開集合にちょうど一致する
——`g`が`P.Ring`の単元であること(`hasMap_X.2`、標準エタール表示の
定義そのもの)から、`g^n`倍は基本開集合を変えない
(`PrimeSpectrum.basicOpen_mul_isUnit`)。`D(f_i)∩D(f_j)`(`Spec B`内)を
`P_i.Ring`の中の基本開集合として実際に特定する、GlueDataの遷移射構成の
核心。

★**sorry 無し**。標準3公理のみ。 -/
theorem standardEtalePresentation_transitionOpen_eq {R S : Type} [CommRing R] [CommRing S]
    [Algebra R S] (Pres : StandardEtalePresentation R S) (x : S) :
    ∃ p : Polynomial R, PrimeSpectrum.basicOpen (Pres.equivRing x) =
      PrimeSpectrum.basicOpen (Polynomial.aeval Pres.P.X p) := by
  obtain ⟨p, n, hpn⟩ := standardEtalePresentation_exists_coord Pres x
  refine ⟨p, ?_⟩
  rw [← hpn, PrimeSpectrum.basicOpen_mul_isUnit _ (Pres.P.hasMap_X.2.pow n)]

end ABC3.Found.CorrHyp
