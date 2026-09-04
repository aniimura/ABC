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

end ABC3.Found.CorrHyp
