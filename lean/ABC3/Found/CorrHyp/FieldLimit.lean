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

open scoped TensorProduct in
/-- **`quotient_mvPolynomial_baseChange`の「純テンソル」上での自然性**
——`(mk I x) ⊗ₜ 1`という(`R`レベルの多項式`x`を単に底変換しただけの)
元の上では、`quotient_mvPolynomial_baseChange`は単に`x`を`MvPolynomial.
map (algebraMap R A)`で送るだけであること。`Algebra.TensorProduct.
quotientTensorEquiv_apply_tmul`(mathlib)+`Ideal.quotientEquiv_mk`
(mathlib)+`quotientMvPolynomialBaseChangeRingEquivAux_comp_algebraMap`
(既存)を繋ぐだけ。`descendPieceR`(`ExtLimit.lean`)を`R`レベルから
さらに昇格した際、局所化のパラメータが正しく対応することを示すための
核心部品——`ideal_map_mvPolynomial_promote_baseChange_eq`と組み合わせて
使う。

★**sorry 無し**。標準3公理のみ。 -/
theorem quotient_mvPolynomial_baseChange_tmul_one (R A ι : Type) [CommRing R] [CommRing A] [Algebra R A]
    (I : Ideal (MvPolynomial ι R)) (x : MvPolynomial ι R) :
    quotient_mvPolynomial_baseChange R A ι I ((Ideal.Quotient.mk I x) ⊗ₜ[R] (1 : A)) =
      Ideal.Quotient.mk (Ideal.map (MvPolynomial.map (algebraMap R A)) I) (MvPolynomial.map (algebraMap R A) x) := by
  unfold quotient_mvPolynomial_baseChange
  simp only [RingEquiv.trans_apply]
  rw [show (Algebra.TensorProduct.quotientTensorEquiv (R := R) R A (MvPolynomial ι R) I).toRingEquiv
      ((Ideal.Quotient.mk I x) ⊗ₜ[R] (1 : A)) =
      (Algebra.TensorProduct.quotientTensorEquiv (R := R) R A (MvPolynomial ι R) I)
        ((Ideal.Quotient.mk I x) ⊗ₜ[R] (1 : A)) from rfl,
    Algebra.TensorProduct.quotientTensorEquiv_apply_tmul, Ideal.quotientEquiv_mk]
  congr 1
  show quotientMvPolynomialBaseChangeRingEquivAux R A ι
      (algebraMap (MvPolynomial ι R) (MvPolynomial ι R ⊗[R] A) x) = MvPolynomial.map (algebraMap R A) x
  exact quotientMvPolynomialBaseChangeRingEquivAux_comp_algebraMap R A ι x

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

open scoped TensorProduct in
/-- **イデアル所属もRレベルへ降りる**——`p`がℝへbase changeした後
`Ideal.span(range map q)`に属するなら、ある共通の精密化`R'`で、
`p`自身(`R'`へ昇格したもの)が`Ideal.span(range 昇格q)`に属する。
`Ideal.mem_span_range_iff_exists_fun`(mathlib、`Fintype`添字の場合の
所属の明示的表示、係数の和)でℝレベルの witness 係数(有限個)を取り出し、
それ自体を`exists_fg_subalgebra_tensor_mvPolynomial_finset`で共通の
`R₀`へ降ろしてから、`exists_mvPolynomial_quotient_ringHom_descend`と
同じ「単射性で等式を`R'`レベルへ押し戻す」パターンを適用する。

`exists_mvPolynomial_quotient_ringHom_descend`が示す「遷移写像の存在の
降下」を、さらに**2つの候補写像の往復が恒等になること**(同型である
ことの確認)まで押し進めるのに使う——遷移**同型**のRレベル降下に必要な
最後のピース。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_mem_ideal_span_range_descend (A : Type) [CommRing A] [Algebra ℚ A]
    (R : FgSubalgebra ℚ ℝ) {ι κ : Type} [Fintype κ]
    (q : κ → MvPolynomial ι (A ⊗[ℚ] R.1)) (p : MvPolynomial ι (A ⊗[ℚ] R.1))
    (h : MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom p ∈
      Ideal.span (Set.range (fun k => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom (q k)))) :
    ∃ (R' : FgSubalgebra ℚ ℝ) (hR : R ≤ R'),
      MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR)).toRingHom p ∈
      Ideal.span (Set.range (fun k => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR)).toRingHom (q k))) := by
  classical
  rw [Ideal.mem_span_range_iff_exists_fun] at h
  obtain ⟨c, hc⟩ := h
  obtain ⟨R₀, hR₀spec⟩ := exists_fg_subalgebra_tensor_mvPolynomial_finset A (Finset.image c Finset.univ)
  choose c₀ hc₀ using fun k : κ => hR₀spec (c k) (Finset.mem_image_of_mem _ (Finset.mem_univ k))
  obtain ⟨R', hR, hR₀⟩ := exists_fgSubalgebra_upperBound2 R R₀
  refine ⟨R', hR, ?_⟩
  rw [Ideal.mem_span_range_iff_exists_fun]
  refine ⟨fun k => MvPolynomial.map
    (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR₀)).toRingHom (c₀ k), ?_⟩
  rw [← sub_eq_zero]
  apply mvPolynomial_algebraTensorMap_val_eq_zero_of_map_eq_zero A R' ι
  rw [map_sub, map_sum, MvPolynomial.map_map, algebraTensorMap_val_comp_inclusion, ← hc]
  have e1 : ∀ k, MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom
      ((MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR₀)).toRingHom (c₀ k)) *
       (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR)).toRingHom (q k))) =
      (c k) * (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom (q k)) := by
    intro k
    rw [map_mul, MvPolynomial.map_map, MvPolynomial.map_map, algebraTensorMap_val_comp_inclusion,
      algebraTensorMap_val_comp_inclusion, hc₀ k]
  rw [Finset.sum_congr rfl (fun k _ => e1 k), sub_self]

open scoped TensorProduct in
/-- **`p`が`Ideal.span(range q)`に属するなら、`R'`へ昇格したものも
`R'`レベルの昇格イデアルに属する**——`Ideal.mem_map_of_mem`(mathlib)+
`Ideal.map_span`を繋ぐだけの単調性。`exists_mem_ideal_span_range_
descend`の逆向き(降下ではなく昇格)版で、複数の`k`ごとに得られる別々の
精密化`R'_k`を共通の`R'`へ揃えるのに使う。

★**sorry 無し**。標準3公理のみ。 -/
theorem mem_ideal_span_range_promote (A : Type) [CommRing A] [Algebra ℚ A]
    {R R' : FgSubalgebra ℚ ℝ} (hR : R ≤ R') {ι κ : Type}
    (q : κ → MvPolynomial ι (A ⊗[ℚ] R.1)) (p : MvPolynomial ι (A ⊗[ℚ] R.1))
    (h : p ∈ Ideal.span (Set.range q)) :
    MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR)).toRingHom p ∈
      Ideal.span (Set.range (fun k => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR)).toRingHom (q k))) := by
  have := Ideal.mem_map_of_mem
    (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR)).toRingHom) h
  rwa [Ideal.map_span, ← Set.range_comp] at this

/-- `Subalgebra.inclusion`を2回昇格させたものは、直接1回で昇格させた
ものと一致すること(`Algebra.TensorProduct.map`で包んだ版)——
`Subalgebra.inclusion_inclusion`(mathlib)+`Algebra.TensorProduct.
map_comp`を繋ぐ。共通の精密化を何段も経由する貼り合わせの配線で
くり返し使う。

★**sorry 無し**。標準3公理のみ。 -/
theorem algebraTensorMap_inclusion_comp_inclusion (A : Type) [CommRing A] [Algebra ℚ A]
    {R S T : FgSubalgebra ℚ ℝ} (h1 : R ≤ S) (h2 : S ≤ T) :
    (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion h2)).toRingHom.comp
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion h1)).toRingHom =
    (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (h1.trans h2))).toRingHom := by
  apply RingHom.ext
  intro x
  show (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion h2))
      ((Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion h1)) x) = _
  rw [← AlgHom.comp_apply, ← Algebra.TensorProduct.map_comp]
  congr 2

/-- `mvPolynomial_map_aeval_comm`の完全に一般な版——係数写像`ψ`が
`Algebra.TensorProduct.map`の形である必要はない(任意の`R→+*S`で
成り立つ)。証明は同一(`MvPolynomial.induction_on`の3ケース)。 -/
theorem mvPolynomial_map_aeval_comm_general {R S : Type} [CommRing R] [CommRing S] (ψ : R →+* S)
    {ι ι' : Type} (p : MvPolynomial ι R) (ev : ι → MvPolynomial ι' R) :
    MvPolynomial.map ψ (MvPolynomial.aeval ev p) =
      MvPolynomial.aeval (fun i => MvPolynomial.map ψ (ev i)) (MvPolynomial.map ψ p) := by
  induction p using MvPolynomial.induction_on with
  | C a =>
    simp only [MvPolynomial.aeval_C, MvPolynomial.algebraMap_eq, MvPolynomial.map_C]
  | add p q hp hq =>
    rw [map_add, map_add, map_add, hp, hq, map_add]
  | mul_X p i hp =>
    simp only [map_mul, MvPolynomial.aeval_X, MvPolynomial.map_X]
    rw [hp]

open scoped TensorProduct in
/-- **`ψ`の値(有限個)だけから、関係式を全く経由せずにRレベルの候補を
取り出す**——`exists_mvPolynomial_quotient_ringHom_descend`から関係式
条件を落としただけの純粋な存在部分。`exists_mvPolynomial_quotient_
ringHom_descend2`で、関係式条件を「等式」ではなく「イデアル所属」へ
弱めるための土台として使う。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_mvPolynomial_eval_descend (A : Type) [CommRing A] [Algebra ℚ A]
    {ι ι' : Type} [Fintype ι] (ψ : ι → MvPolynomial ι' (A ⊗[ℚ] ℝ)) :
    ∃ (R' : FgSubalgebra ℚ ℝ) (ev : ι → MvPolynomial ι' (A ⊗[ℚ] R'.1)),
      ∀ i, MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom (ev i)
        = ψ i := by
  classical
  obtain ⟨R', hR'spec⟩ := exists_fg_subalgebra_tensor_mvPolynomial_finset A (Finset.image ψ Finset.univ)
  choose ev hev using fun i : ι => hR'spec (ψ i) (Finset.mem_image_of_mem _ (Finset.mem_univ i))
  exact ⟨R', ev, hev⟩

open scoped TensorProduct in
/-- **「項目(d)の第二段」の完成形——2つの異なるRレベル片(異なる`R`・
`R₂`)の間の遷移写像が、共通の精密化`R'`へ構成的に降りる**。
`exists_mvPolynomial_quotient_ringHom_descend`の一般化——関係式条件を
「等式`=0`」ではなく「イデアル所属」で述べることで、**目標側も
quotient(genuine な片同士)である場合**に対応した。

`ψ`(既知のℝレベルの遷移写像が生成元に送る値)の存在部分
(`exists_mvPolynomial_eval_descend`)と関係式の所属部分
(`mvPolynomial_map_aeval_comm_general`+`exists_mem_ideal_span_range_
descend`)を独立に降ろしてから、有限個の`k`ごとに異なりうる精密化を
`exists_fgSubalgebra_upperBound`で単一の`R'`へ合流させる
(`mem_ideal_span_range_promote`で個別の結果を`R'`へ昇格)。

これで、2つの隣接する`descendPieceR`(異なる`R_i`)の重なり部分の
遷移写像が、実際に共通のRレベル精密化上の候補写像として構成的に
得られることが証明された——`transitionElem`/`gdT`/`cocycle`のRレベル
版に相当する遷移データの降下の核心が完成した。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_mvPolynomial_quotient_ringHom_descend2 (A : Type) [CommRing A] [Algebra ℚ A]
    (R R₂ : FgSubalgebra ℚ ℝ) {ι κ ι' κ' : Type} [Fintype ι] [Fintype κ] [Fintype κ']
    (q : κ → MvPolynomial ι (A ⊗[ℚ] R.1)) (q₂ : κ' → MvPolynomial ι' (A ⊗[ℚ] R₂.1))
    (ψ : ι → MvPolynomial ι' (A ⊗[ℚ] ℝ))
    (hψ : ∀ k, MvPolynomial.aeval ψ
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom (q k)) ∈
      Ideal.span (Set.range (fun k' => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R₂.1)).toRingHom (q₂ k')))) :
    ∃ (R' : FgSubalgebra ℚ ℝ) (hR : R ≤ R') (hR₂ : R₂ ≤ R') (ev : ι → MvPolynomial ι' (A ⊗[ℚ] R'.1)),
      (∀ i, MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom (ev i)
        = ψ i) ∧
      (∀ k, MvPolynomial.aeval ev (MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR)).toRingHom (q k)) ∈
        Ideal.span (Set.range (fun k' => MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR₂)).toRingHom (q₂ k')))) := by
  classical
  obtain ⟨R₀, ev₀, hev₀⟩ := exists_mvPolynomial_eval_descend A ψ
  obtain ⟨Ra, hRRa, hR₂Ra⟩ := exists_fgSubalgebra_upperBound2 R R₂
  obtain ⟨R₁, hRaR₁, hR₀R₁⟩ := exists_fgSubalgebra_upperBound2 Ra R₀
  have hR : R ≤ R₁ := hRRa.trans hRaR₁
  have hR₂ : R₂ ≤ R₁ := hR₂Ra.trans hRaR₁
  obtain ⟨ev₁, hev₁⟩ : ∃ ev₁ : ι → MvPolynomial ι' (A ⊗[ℚ] R₁.1), ∀ i, MvPolynomial.map
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R₁.1)).toRingHom (ev₁ i) = ψ i :=
    ⟨fun i => MvPolynomial.map
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR₀R₁)).toRingHom (ev₀ i), fun i => by
      rw [MvPolynomial.map_map, algebraTensorMap_val_comp_inclusion]; exact hev₀ i⟩
  have hcombine : ∀ k, MvPolynomial.map
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R₁.1)).toRingHom
      (MvPolynomial.aeval ev₁ (MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR)).toRingHom (q k))) =
      MvPolynomial.aeval ψ
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom (q k)) := by
    intro k
    rw [mvPolynomial_map_aeval_comm]
    rw [show (fun i => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R₁.1)).toRingHom (ev₁ i)) = ψ from
      funext hev₁]
    rw [MvPolynomial.map_map, algebraTensorMap_val_comp_inclusion]
  have hmemR₁ : ∀ k, MvPolynomial.map
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R₁.1)).toRingHom
      (MvPolynomial.aeval ev₁ (MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR)).toRingHom (q k))) ∈
      Ideal.span (Set.range (fun k' => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R₁.1)).toRingHom
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR₂)).toRingHom
          (q₂ k')))) := by
    intro k
    rw [hcombine]
    have hspan : Set.range (fun k' => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R₁.1)).toRingHom
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR₂)).toRingHom
          (q₂ k'))) = Set.range (fun k' => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R₂.1)).toRingHom (q₂ k')) := by
      apply congrArg Set.range
      funext k'
      rw [MvPolynomial.map_map, algebraTensorMap_val_comp_inclusion]
    rw [hspan]
    exact hψ k
  choose Rk hR₁Rk hmemk using fun k => exists_mem_ideal_span_range_descend A R₁
    (fun k' : κ' => MvPolynomial.map
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR₂)).toRingHom (q₂ k'))
    (MvPolynomial.aeval ev₁ (MvPolynomial.map
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR)).toRingHom (q k)))
    (hmemR₁ k)
  obtain ⟨R₀', hR₀'spec⟩ := exists_fgSubalgebra_upperBound (Finset.univ : Finset κ) Rk
  obtain ⟨R', hR₁R', hR₀'R'⟩ := exists_fgSubalgebra_upperBound2 R₁ R₀'
  have hRkR' : ∀ k, Rk k ≤ R' := fun k => (hR₀'spec k (Finset.mem_univ k)).trans hR₀'R'
  refine ⟨R', hR.trans hR₁R', hR₂.trans hR₁R', fun i => MvPolynomial.map
    (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR₁R')).toRingHom (ev₁ i), ?_, ?_⟩
  · intro i
    rw [MvPolynomial.map_map, algebraTensorMap_val_comp_inclusion]
    exact hev₁ i
  · intro k
    have hpromote := mem_ideal_span_range_promote A (hRkR' k) _ _ (hmemk k)
    rw [Ideal.mem_span_range_iff_exists_fun] at hpromote ⊢
    obtain ⟨c, hc⟩ := hpromote
    refine ⟨c, ?_⟩
    have e1 : ∀ k', MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hR₂.trans hR₁R'))).toRingHom (q₂ k') =
        MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hRkR' k))).toRingHom
          (MvPolynomial.map
            (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hR₁Rk k))).toRingHom
            (MvPolynomial.map
              (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR₂)).toRingHom (q₂ k'))) := by
      intro k'
      rw [MvPolynomial.map_map, MvPolynomial.map_map, algebraTensorMap_inclusion_comp_inclusion,
        algebraTensorMap_inclusion_comp_inclusion]
    have e2 : MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hR.trans hR₁R'))).toRingHom (q k) =
        MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hRkR' k))).toRingHom
          (MvPolynomial.map
            (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hR₁Rk k))).toRingHom
            (MvPolynomial.map
              (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR)).toRingHom (q k))) := by
      rw [MvPolynomial.map_map, MvPolynomial.map_map, algebraTensorMap_inclusion_comp_inclusion,
        algebraTensorMap_inclusion_comp_inclusion]
    have e3 : MvPolynomial.aeval (fun i => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR₁R')).toRingHom (ev₁ i))
        (MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hR.trans hR₁R'))).toRingHom (q k)) =
        MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hRkR' k))).toRingHom
          (MvPolynomial.map
            (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hR₁Rk k))).toRingHom
            (MvPolynomial.aeval ev₁ (MvPolynomial.map
              (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR)).toRingHom (q k)))) := by
      rw [e2, mvPolynomial_map_aeval_comm_general, mvPolynomial_map_aeval_comm_general]
      have hevfun : (fun i => MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR₁R')).toRingHom (ev₁ i)) =
        (fun i => MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hRkR' k))).toRingHom
          (MvPolynomial.map
            (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hR₁Rk k))).toRingHom (ev₁ i))) := by
        funext i
        rw [MvPolynomial.map_map, algebraTensorMap_inclusion_comp_inclusion]
      rw [hevfun]
    simp only [e1]
    rw [e3, hc]

open scoped TensorProduct in
/-- `Subalgebra.inclusion`を2回昇格させた後で生成元代入をとったものと、
1回で昇格させてから生成元代入をとったものが一致すること
(`round_trip_promote_eq2`の非対称版・準備補題)——遷移同型の往復合成
`aeval ev'(ev i)`をさらに1段昇格させる場面で使う。

★**sorry 無し**。標準3公理のみ。 -/
theorem round_trip_promote_eq (A : Type) [CommRing A] [Algebra ℚ A]
    (Rc Rk : FgSubalgebra ℚ ℝ) (hRck : Rc ≤ Rk) (R' : FgSubalgebra ℚ ℝ) (hRkR' : Rk ≤ R') (hRcR' : Rc ≤ R')
    {ι ι' : Type} (ev : ι → MvPolynomial ι' (A ⊗[ℚ] Rc.1)) (ev' : ι' → MvPolynomial ι (A ⊗[ℚ] Rc.1)) (i : ι) :
    MvPolynomial.aeval (fun j => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRcR')).toRingHom (ev' j))
      (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRcR')).toRingHom
        (ev i)) - MvPolynomial.X i =
      MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRkR')).toRingHom
        (MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRck)).toRingHom
          (MvPolynomial.aeval ev' (ev i) - MvPolynomial.X i)) := by
  have step1 : MvPolynomial.map
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRkR')).toRingHom
      (MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRck)).toRingHom
        (MvPolynomial.aeval ev' (ev i))) =
      MvPolynomial.aeval (fun j => MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRcR')).toRingHom (ev' j))
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRcR')).toRingHom
          (ev i)) := by
    rw [MvPolynomial.map_map, algebraTensorMap_inclusion_comp_inclusion, mvPolynomial_map_aeval_comm_general]
  have step2 : MvPolynomial.map
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRkR')).toRingHom
      (MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRck)).toRingHom
        (MvPolynomial.X (σ := ι) i)) = MvPolynomial.X i := by
    rw [MvPolynomial.map_map, algebraTensorMap_inclusion_comp_inclusion, MvPolynomial.map_X]
  rw [map_sub, map_sub, step1, step2]

open scoped TensorProduct in
/-- **候補となる遷移写像の往復合成が、途中の精密化を経由しても等しい
こと**——`ev`・`ev'`(共通の`Rc`レベル)を`Rck:Rc→Rk`経由でさらに昇格
させた場合と、直接`Rc→R'`(`hRcR'`)で昇格させた場合とで、往復合成
`aeval ev'(ev i) - X i`が一致する。`exists_round_trip_descend`が返す
結果を、さらに別の精密化`R'`へ持ち上げる際に使う——「Rレベルの
候補写像の往復が恒等であること」を最終的な共通段階`R'`へ集約する
配線の核心。

★**sorry 無し**。標準3公理のみ。 -/
theorem round_trip_promote_eq2 (A : Type) [CommRing A] [Algebra ℚ A]
    (Rc Rk : FgSubalgebra ℚ ℝ) (hRck : Rc ≤ Rk) (R' : FgSubalgebra ℚ ℝ) (hRkR' : Rk ≤ R') (hRcR' : Rc ≤ R')
    {ι ι' : Type} (ev : ι → MvPolynomial ι' (A ⊗[ℚ] Rc.1)) (ev' : ι' → MvPolynomial ι (A ⊗[ℚ] Rc.1)) (i : ι) :
    MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRkR')).toRingHom
      (MvPolynomial.aeval (fun j => MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRck)).toRingHom (ev' j))
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRck)).toRingHom
          (ev i)) - MvPolynomial.X i) =
      MvPolynomial.aeval (fun j => MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRcR')).toRingHom (ev' j))
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRcR')).toRingHom
          (ev i)) - MvPolynomial.X i := by
  have step1 : MvPolynomial.map
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRkR')).toRingHom
      (MvPolynomial.aeval (fun j => MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRck)).toRingHom (ev' j))
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRck)).toRingHom
          (ev i))) =
      MvPolynomial.aeval (fun j => MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRcR')).toRingHom (ev' j))
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRcR')).toRingHom
          (ev i)) := by
    rw [mvPolynomial_map_aeval_comm_general, MvPolynomial.map_map, algebraTensorMap_inclusion_comp_inclusion]
    congr 1
    congr 1
    funext j
    rw [MvPolynomial.map_map, algebraTensorMap_inclusion_comp_inclusion]
  rw [map_sub, MvPolynomial.map_X, step1]

open scoped TensorProduct in
/-- **遷移写像の候補の往復合成が、既知のℝレベルの恒等性からRレベルへ
構成的に降りる**——`ψ`・`ψ'`(ℝレベル)の往復合成`aeval ψ'(ψ i) - X i`が
`q`の関係式イデアルに属するなら(既知)、`ev`・`ev'`(共通の`Rc`レベル、
`ψ`・`ψ'`をそれぞれ再現する)自身の往復合成も、ある共通のさらなる
精密化`R'`で同じイデアルに属する。`exists_mem_ideal_span_range_
descend`(有限個の`i`ごとに個別に適用)+`exists_fgSubalgebra_
upperBound`(それらの精密化を単一の`R'`へ合流)を組み合わせる。

`exists_mvPolynomial_quotient_ringEquiv_descend`で、遷移写像の候補が
実際に同型であることを示す核心部品の1つ。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_round_trip_descend (A : Type) [CommRing A] [Algebra ℚ A]
    (R Rc : FgSubalgebra ℚ ℝ) (hR : R ≤ Rc) {ι ι' κ : Type} [Fintype ι] [Fintype κ]
    (q : κ → MvPolynomial ι (A ⊗[ℚ] R.1))
    (ev : ι → MvPolynomial ι' (A ⊗[ℚ] Rc.1)) (ev' : ι' → MvPolynomial ι (A ⊗[ℚ] Rc.1))
    (ψ : ι → MvPolynomial ι' (A ⊗[ℚ] ℝ)) (ψ' : ι' → MvPolynomial ι (A ⊗[ℚ] ℝ))
    (hev : ∀ i, MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val Rc.1)).toRingHom
        (ev i) = ψ i)
    (hev' : ∀ j, MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val Rc.1)).toRingHom
        (ev' j) = ψ' j)
    (hround : ∀ i, MvPolynomial.aeval ψ' (ψ i) - MvPolynomial.X i ∈
      Ideal.span (Set.range (fun k => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom (q k)))) :
    ∃ (R' : FgSubalgebra ℚ ℝ) (hRc : Rc ≤ R'),
      ∀ i, MvPolynomial.aeval (fun j => MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRc)).toRingHom (ev' j))
          (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRc)).toRingHom
            (ev i)) - MvPolynomial.X i ∈
        Ideal.span (Set.range (fun k => MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hR.trans hRc))).toRingHom (q k))) := by
  classical
  have hmemRc : ∀ i, MvPolynomial.map
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val Rc.1)).toRingHom
      (MvPolynomial.aeval ev' (ev i) - MvPolynomial.X i) ∈
      Ideal.span (Set.range (fun k => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom (q k))) := by
    intro i
    rw [map_sub, mvPolynomial_map_aeval_comm]
    rw [show (fun j => MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val Rc.1)).toRingHom
        (ev' j)) = ψ' from funext hev', hev, MvPolynomial.map_X]
    exact hround i
  have hmemRc' : ∀ i, MvPolynomial.map
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val Rc.1)).toRingHom
      (MvPolynomial.aeval ev' (ev i) - MvPolynomial.X i) ∈
      Ideal.span (Set.range (fun k => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val Rc.1)).toRingHom
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR)).toRingHom
          (q k)))) := by
    intro i
    have hspan : Set.range (fun k => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val Rc.1)).toRingHom
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR)).toRingHom
          (q k))) = Set.range (fun k => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom (q k)) := by
      apply congrArg Set.range
      funext k
      rw [MvPolynomial.map_map, algebraTensorMap_val_comp_inclusion]
    rw [hspan]
    exact hmemRc i
  choose Rk hRck hmemk using fun i => exists_mem_ideal_span_range_descend A Rc
    (fun k => MvPolynomial.map
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR)).toRingHom (q k))
    (MvPolynomial.aeval ev' (ev i) - MvPolynomial.X i) (hmemRc' i)
  obtain ⟨R₀', hR₀'spec⟩ := exists_fgSubalgebra_upperBound (Finset.univ : Finset ι) Rk
  obtain ⟨R', hRcR', hR₀'R'⟩ := exists_fgSubalgebra_upperBound2 Rc R₀'
  have hRkR' : ∀ i, Rk i ≤ R' := fun i => (hR₀'spec i (Finset.mem_univ i)).trans hR₀'R'
  refine ⟨R', hRcR', ?_⟩
  intro i
  have hpromote := mem_ideal_span_range_promote A (hRkR' i) _ _ (hmemk i)
  rw [Ideal.mem_span_range_iff_exists_fun] at hpromote ⊢
  obtain ⟨c, hc⟩ := hpromote
  refine ⟨c, ?_⟩
  have e1 : ∀ k, MvPolynomial.map
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hR.trans hRcR'))).toRingHom (q k) =
      MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hRkR' i))).toRingHom
        (MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hRck i))).toRingHom
          (MvPolynomial.map
            (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR)).toRingHom (q k))) := by
    intro k
    rw [MvPolynomial.map_map, MvPolynomial.map_map, algebraTensorMap_inclusion_comp_inclusion,
      algebraTensorMap_inclusion_comp_inclusion]
  simp only [e1]
  rw [round_trip_promote_eq A Rc (Rk i) (hRck i) R' (hRkR' i) hRcR' ev ev' i, hc]

set_option maxHeartbeats 4000000 in
open scoped TensorProduct in
/-- **「候補の遷移写像が実際に同型であること」の完成——項目(d)の
第二段の最終形**: `R`レベルの片1(`q`、生成元`ι`)と`R₂`レベルの片2
(`q₂`、生成元`ι'`)の間の、既知のℝレベルの相互に逆な遷移写像
(`ψ`・`ψ'`、往復合成が恒等になる)から出発し、ある共通の精密化`R'`上の
候補写像`ev`・`ev'`が、`ψ`・`ψ'`を再現し、互いの関係式を互いのイデアル
へ写し、かつ**往復合成が実際に恒等射であること**まで構成的に示す。

`exists_mvPolynomial_quotient_ringHom_descend2`を両方向(`ψ`用・`ψ'`用)
に適用して`ev`・`ev'`の存在部分と関係式部分を得たのち、共通の精密化
`Rc`へ合流させ、その上で`exists_round_trip_descend`を両方向(`hround1`・
`hround2`用)に適用して往復合成の恒等性を得て、最後にすべてを単一の
`R'`へ揃える——`transitionElem`/`gdT`/`cocycle`のRレベル版に相当する
遷移**同型**データの降下が、これで完全に構成的に証明された。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_mvPolynomial_quotient_ringEquiv_descend (A : Type) [CommRing A] [Algebra ℚ A]
    (R R₂ : FgSubalgebra ℚ ℝ) {ι κ ι' κ' : Type} [Fintype ι] [Fintype ι'] [Fintype κ] [Fintype κ']
    (q : κ → MvPolynomial ι (A ⊗[ℚ] R.1)) (q₂ : κ' → MvPolynomial ι' (A ⊗[ℚ] R₂.1))
    (ψ : ι → MvPolynomial ι' (A ⊗[ℚ] ℝ)) (ψ' : ι' → MvPolynomial ι (A ⊗[ℚ] ℝ))
    (hψ : ∀ k, MvPolynomial.aeval ψ
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom (q k)) ∈
      Ideal.span (Set.range (fun k' => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R₂.1)).toRingHom (q₂ k'))))
    (hψ' : ∀ k', MvPolynomial.aeval ψ'
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R₂.1)).toRingHom (q₂ k')) ∈
      Ideal.span (Set.range (fun k => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom (q k))))
    (hround1 : ∀ i, MvPolynomial.aeval ψ' (ψ i) - MvPolynomial.X i ∈
      Ideal.span (Set.range (fun k => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom (q k))))
    (hround2 : ∀ j, MvPolynomial.aeval ψ (ψ' j) - MvPolynomial.X j ∈
      Ideal.span (Set.range (fun k' => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R₂.1)).toRingHom (q₂ k')))) :
    ∃ (R' : FgSubalgebra ℚ ℝ) (hR : R ≤ R') (hR₂ : R₂ ≤ R')
      (ev : ι → MvPolynomial ι' (A ⊗[ℚ] R'.1)) (ev' : ι' → MvPolynomial ι (A ⊗[ℚ] R'.1)),
      (∀ i, MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom (ev i)
        = ψ i) ∧
      (∀ j, MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom (ev' j)
        = ψ' j) ∧
      (∀ k, MvPolynomial.aeval ev (MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR)).toRingHom (q k)) ∈
        Ideal.span (Set.range (fun k' => MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR₂)).toRingHom (q₂ k')))) ∧
      (∀ k', MvPolynomial.aeval ev' (MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR₂)).toRingHom (q₂ k')) ∈
        Ideal.span (Set.range (fun k => MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR)).toRingHom (q k)))) ∧
      (∀ i, MvPolynomial.aeval ev' (ev i) - MvPolynomial.X i ∈
        Ideal.span (Set.range (fun k => MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR)).toRingHom (q k)))) ∧
      (∀ j, MvPolynomial.aeval ev (ev' j) - MvPolynomial.X j ∈
        Ideal.span (Set.range (fun k' => MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR₂)).toRingHom (q₂ k')))) := by
  classical
  obtain ⟨Ra, hRRa, hR₂Ra, ev_a, hev_a, hq_a⟩ := exists_mvPolynomial_quotient_ringHom_descend2 A R R₂ q q₂ ψ hψ
  obtain ⟨Rb, hR₂Rb, hRRb, ev'_b, hev'_b, hq₂_b⟩ := exists_mvPolynomial_quotient_ringHom_descend2 A R₂ R q₂ q ψ' hψ'
  obtain ⟨Rc, hRaRc, hRbRc⟩ := exists_fgSubalgebra_upperBound2 Ra Rb
  set ev_c : ι → MvPolynomial ι' (A ⊗[ℚ] Rc.1) := fun i => MvPolynomial.map
    (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRaRc)).toRingHom (ev_a i) with hev_c_def
  set ev'_c : ι' → MvPolynomial ι (A ⊗[ℚ] Rc.1) := fun j => MvPolynomial.map
    (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRbRc)).toRingHom (ev'_b j) with hev'_c_def
  have hev_c : ∀ i, MvPolynomial.map
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val Rc.1)).toRingHom (ev_c i) = ψ i := by
    intro i
    rw [hev_c_def, MvPolynomial.map_map, algebraTensorMap_val_comp_inclusion]
    exact hev_a i
  have hev'_c : ∀ j, MvPolynomial.map
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val Rc.1)).toRingHom (ev'_c j) = ψ' j := by
    intro j
    rw [hev'_c_def, MvPolynomial.map_map, algebraTensorMap_val_comp_inclusion]
    exact hev'_b j
  have hq_c : ∀ k, MvPolynomial.aeval ev_c (MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hRRa.trans hRaRc))).toRingHom (q k)) ∈
      Ideal.span (Set.range (fun k' => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hR₂Ra.trans hRaRc))).toRingHom
        (q₂ k'))) := by
    intro k
    have hp := mem_ideal_span_range_promote A hRaRc _ _ (hq_a k)
    have eL : MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRaRc)).toRingHom
        (MvPolynomial.aeval ev_a (MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRRa)).toRingHom (q k))) =
        MvPolynomial.aeval ev_c (MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hRRa.trans hRaRc))).toRingHom (q k)) := by
      rw [mvPolynomial_map_aeval_comm_general, MvPolynomial.map_map, algebraTensorMap_inclusion_comp_inclusion]
    have eR : Set.range (fun k' => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRaRc)).toRingHom
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR₂Ra)).toRingHom
          (q₂ k'))) = Set.range (fun k' => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hR₂Ra.trans hRaRc))).toRingHom
        (q₂ k')) := by
      apply congrArg Set.range
      funext k'
      rw [MvPolynomial.map_map, algebraTensorMap_inclusion_comp_inclusion]
    rw [eL, eR] at hp
    exact hp
  have hq₂_c : ∀ k', MvPolynomial.aeval ev'_c (MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hR₂Rb.trans hRbRc))).toRingHom
        (q₂ k')) ∈
      Ideal.span (Set.range (fun k => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hRRb.trans hRbRc))).toRingHom
        (q k))) := by
    intro k'
    have hp := mem_ideal_span_range_promote A hRbRc _ _ (hq₂_b k')
    have eL : MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRbRc)).toRingHom
        (MvPolynomial.aeval ev'_b (MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR₂Rb)).toRingHom (q₂ k'))) =
        MvPolynomial.aeval ev'_c (MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hR₂Rb.trans hRbRc))).toRingHom
          (q₂ k')) := by
      rw [mvPolynomial_map_aeval_comm_general, MvPolynomial.map_map, algebraTensorMap_inclusion_comp_inclusion]
    have eR : Set.range (fun k => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRbRc)).toRingHom
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hRRb)).toRingHom
          (q k))) = Set.range (fun k => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hRRb.trans hRbRc))).toRingHom
        (q k)) := by
      apply congrArg Set.range
      funext k
      rw [MvPolynomial.map_map, algebraTensorMap_inclusion_comp_inclusion]
    rw [eL, eR] at hp
    exact hp
  obtain ⟨R'₁, hRcR'₁, hround1'⟩ := exists_round_trip_descend A R Rc (hRRa.trans hRaRc) q ev_c ev'_c ψ ψ' hev_c hev'_c
    hround1
  obtain ⟨R'₂, hRcR'₂, hround2'⟩ := exists_round_trip_descend A R₂ Rc (hR₂Rb.trans hRbRc) q₂ ev'_c ev_c ψ' ψ hev'_c
    hev_c hround2
  obtain ⟨R', hR'₁R', hR'₂R'⟩ := exists_fgSubalgebra_upperBound2 R'₁ R'₂
  refine ⟨R', hRRa.trans (hRaRc.trans (hRcR'₁.trans hR'₁R')), hR₂Rb.trans (hRbRc.trans (hRcR'₂.trans hR'₂R')),
    fun i => MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A)
      (Subalgebra.inclusion (hRcR'₁.trans hR'₁R'))).toRingHom (ev_c i),
    fun j => MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A)
      (Subalgebra.inclusion (hRcR'₁.trans hR'₁R'))).toRingHom (ev'_c j), ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro i
    rw [MvPolynomial.map_map, algebraTensorMap_val_comp_inclusion]
    exact hev_c i
  · intro j
    rw [MvPolynomial.map_map, algebraTensorMap_val_comp_inclusion]
    exact hev'_c j
  · intro k
    have hp := mem_ideal_span_range_promote A (hRcR'₁.trans hR'₁R') _ _ (hq_c k)
    have eL : MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A)
          (Subalgebra.inclusion (hRcR'₁.trans hR'₁R'))).toRingHom
        (MvPolynomial.aeval ev_c (MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hRRa.trans hRaRc))).toRingHom (q k))) =
        MvPolynomial.aeval (fun i => MvPolynomial.map
            (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hRcR'₁.trans hR'₁R'))).toRingHom
            (ev_c i))
          (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A)
            (Subalgebra.inclusion (hRRa.trans (hRaRc.trans (hRcR'₁.trans hR'₁R'))))).toRingHom (q k)) := by
      rw [mvPolynomial_map_aeval_comm_general, MvPolynomial.map_map, algebraTensorMap_inclusion_comp_inclusion]
    have eR : Set.range (fun k' => MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A)
          (Subalgebra.inclusion (hRcR'₁.trans hR'₁R'))).toRingHom
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A)
          (Subalgebra.inclusion (hR₂Ra.trans hRaRc))).toRingHom (q₂ k'))) = Set.range (fun k' =>
        MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A)
          (Subalgebra.inclusion (hR₂Rb.trans (hRbRc.trans (hRcR'₂.trans hR'₂R'))))).toRingHom (q₂ k')) := by
      apply congrArg Set.range
      funext k'
      rw [MvPolynomial.map_map, algebraTensorMap_inclusion_comp_inclusion]
    rw [eL, eR] at hp
    exact hp
  · intro k'
    have hp := mem_ideal_span_range_promote A (hRcR'₂.trans hR'₂R') _ _ (hq₂_c k')
    have eL : MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A)
          (Subalgebra.inclusion (hRcR'₂.trans hR'₂R'))).toRingHom
        (MvPolynomial.aeval ev'_c (MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hR₂Rb.trans hRbRc))).toRingHom
          (q₂ k'))) =
        MvPolynomial.aeval (fun j => MvPolynomial.map
            (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion (hRcR'₁.trans hR'₁R'))).toRingHom
            (ev'_c j))
          (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A)
            (Subalgebra.inclusion (hR₂Rb.trans (hRbRc.trans (hRcR'₂.trans hR'₂R'))))).toRingHom (q₂ k')) := by
      rw [mvPolynomial_map_aeval_comm_general, MvPolynomial.map_map, algebraTensorMap_inclusion_comp_inclusion]
    have eR : Set.range (fun k => MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A)
          (Subalgebra.inclusion (hRcR'₂.trans hR'₂R'))).toRingHom
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A)
          (Subalgebra.inclusion (hRRb.trans hRbRc))).toRingHom (q k))) = Set.range (fun k =>
        MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A)
          (Subalgebra.inclusion (hRRa.trans (hRaRc.trans (hRcR'₁.trans hR'₁R'))))).toRingHom (q k)) := by
      apply congrArg Set.range
      funext k
      rw [MvPolynomial.map_map, algebraTensorMap_inclusion_comp_inclusion]
    rw [eL, eR] at hp
    exact hp
  · intro i
    have hp := mem_ideal_span_range_promote A hR'₁R' _ _ (hround1' i)
    rw [round_trip_promote_eq2 A Rc R'₁ hRcR'₁ R' hR'₁R' (hRcR'₁.trans hR'₁R') ev_c ev'_c i] at hp
    have eR : Set.range (fun k => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR'₁R')).toRingHom
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A)
          (Subalgebra.inclusion ((hRRa.trans hRaRc).trans hRcR'₁))).toRingHom (q k))) = Set.range (fun k =>
        MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A)
          (Subalgebra.inclusion (hRRa.trans (hRaRc.trans (hRcR'₁.trans hR'₁R'))))).toRingHom (q k)) := by
      apply congrArg Set.range
      funext k
      rw [MvPolynomial.map_map, algebraTensorMap_inclusion_comp_inclusion]
    rwa [eR] at hp
  · intro j
    have hp := mem_ideal_span_range_promote A hR'₂R' _ _ (hround2' j)
    rw [round_trip_promote_eq2 A Rc R'₂ hRcR'₂ R' hR'₂R' (hRcR'₂.trans hR'₂R') ev'_c ev_c j] at hp
    have eR : Set.range (fun k' => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR'₂R')).toRingHom
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A)
          (Subalgebra.inclusion ((hR₂Rb.trans hRbRc).trans hRcR'₂))).toRingHom (q₂ k'))) = Set.range (fun k' =>
        MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A)
          (Subalgebra.inclusion (hR₂Rb.trans (hRbRc.trans (hRcR'₂.trans hR'₂R'))))).toRingHom (q₂ k')) := by
      apply congrArg Set.range
      funext k'
      rw [MvPolynomial.map_map, algebraTensorMap_inclusion_comp_inclusion]
    rwa [eR] at hp

/-- **`exists_mvPolynomial_quotient_ringEquiv_descend`の生データから
実際の`RingEquiv`を組み立てる**——`Ideal.Quotient.lift`で`ev`・`ev'`
それぞれから商環の間のRingHom `f`・`g`を構成し(well-definedness は
`hq`・`hq₂`から`Ideal.span_le`経由)、`RingEquiv.ofRingHom`で束ねる。
往復が恒等であること(`f.comp g = id`・`g.comp f = id`)は
`Ideal.Quotient.ringHom_ext`(商からの一致は`mk`との合成で十分)+
`MvPolynomial.ringHom_ext`(`C`・`X`上の一致で十分)で`X`の場合に
`hround1`・`hround2`(`Ideal.Quotient.eq`経由)を、`C`の場合に
`aeval_C`(定数を保つ)を使うだけで閉じる。汎用的な純代数の補題
——CorrHyp固有のデータには依存しない。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_mvPolynomial_quotient_ringEquiv_of_data {T : Type} [CommRing T]
    {ι κ ι' κ' : Type} (q : κ → MvPolynomial ι T) (q₂ : κ' → MvPolynomial ι' T)
    (ev : ι → MvPolynomial ι' T) (ev' : ι' → MvPolynomial ι T)
    (hq : ∀ k, MvPolynomial.aeval ev (q k) ∈ Ideal.span (Set.range q₂))
    (hq₂ : ∀ k', MvPolynomial.aeval ev' (q₂ k') ∈ Ideal.span (Set.range q))
    (hround1 : ∀ i, MvPolynomial.aeval ev' (ev i) - MvPolynomial.X i ∈ Ideal.span (Set.range q))
    (hround2 : ∀ j, MvPolynomial.aeval ev (ev' j) - MvPolynomial.X j ∈ Ideal.span (Set.range q₂)) :
    Nonempty ((MvPolynomial ι T ⧸ Ideal.span (Set.range q)) ≃+*
      (MvPolynomial ι' T ⧸ Ideal.span (Set.range q₂))) := by
  have hwd1 : ∀ a ∈ Ideal.span (Set.range q), (Ideal.Quotient.mk (Ideal.span (Set.range q₂)))
      ((MvPolynomial.aeval ev).toRingHom a) = 0 := by
    intro a ha
    rw [Ideal.Quotient.eq_zero_iff_mem]
    have hle : Ideal.span (Set.range q) ≤
        Ideal.comap (MvPolynomial.aeval ev).toRingHom (Ideal.span (Set.range q₂)) :=
      Ideal.span_le.mpr (fun x hx => by obtain ⟨k, rfl⟩ := hx; exact hq k)
    exact hle ha
  have hwd2 : ∀ a ∈ Ideal.span (Set.range q₂), (Ideal.Quotient.mk (Ideal.span (Set.range q)))
      ((MvPolynomial.aeval ev').toRingHom a) = 0 := by
    intro a ha
    rw [Ideal.Quotient.eq_zero_iff_mem]
    have hle : Ideal.span (Set.range q₂) ≤
        Ideal.comap (MvPolynomial.aeval ev').toRingHom (Ideal.span (Set.range q)) :=
      Ideal.span_le.mpr (fun x hx => by obtain ⟨k, rfl⟩ := hx; exact hq₂ k)
    exact hle ha
  set f := Ideal.Quotient.lift (Ideal.span (Set.range q))
    ((Ideal.Quotient.mk (Ideal.span (Set.range q₂))).comp (MvPolynomial.aeval ev).toRingHom) hwd1 with hf_def
  set g := Ideal.Quotient.lift (Ideal.span (Set.range q₂))
    ((Ideal.Quotient.mk (Ideal.span (Set.range q))).comp (MvPolynomial.aeval ev').toRingHom) hwd2 with hg_def
  refine ⟨RingEquiv.ofRingHom f g ?_ ?_⟩
  · apply Ideal.Quotient.ringHom_ext
    apply MvPolynomial.ringHom_ext
    · intro a
      show f (g (Ideal.Quotient.mk (Ideal.span (Set.range q₂)) (MvPolynomial.C a))) =
        Ideal.Quotient.mk (Ideal.span (Set.range q₂)) (MvPolynomial.C a)
      rw [hg_def, Ideal.Quotient.lift_mk]
      show f (Ideal.Quotient.mk (Ideal.span (Set.range q)) (MvPolynomial.aeval ev' (MvPolynomial.C a))) = _
      rw [MvPolynomial.aeval_C, hf_def, Ideal.Quotient.lift_mk]
      show Ideal.Quotient.mk (Ideal.span (Set.range q₂)) (MvPolynomial.aeval ev (MvPolynomial.C a)) = _
      rw [MvPolynomial.aeval_C]
      rfl
    · intro i
      show f (g (Ideal.Quotient.mk (Ideal.span (Set.range q₂)) (MvPolynomial.X i))) =
        Ideal.Quotient.mk (Ideal.span (Set.range q₂)) (MvPolynomial.X i)
      rw [hg_def, Ideal.Quotient.lift_mk]
      show f (Ideal.Quotient.mk (Ideal.span (Set.range q)) (MvPolynomial.aeval ev' (MvPolynomial.X i))) = _
      rw [MvPolynomial.aeval_X, hf_def, Ideal.Quotient.lift_mk]
      show Ideal.Quotient.mk (Ideal.span (Set.range q₂)) (MvPolynomial.aeval ev (ev' i)) = _
      rw [Ideal.Quotient.eq]
      exact hround2 i
  · apply Ideal.Quotient.ringHom_ext
    apply MvPolynomial.ringHom_ext
    · intro a
      show g (f (Ideal.Quotient.mk (Ideal.span (Set.range q)) (MvPolynomial.C a))) =
        Ideal.Quotient.mk (Ideal.span (Set.range q)) (MvPolynomial.C a)
      rw [hf_def, Ideal.Quotient.lift_mk]
      show g (Ideal.Quotient.mk (Ideal.span (Set.range q₂)) (MvPolynomial.aeval ev (MvPolynomial.C a))) = _
      rw [MvPolynomial.aeval_C, hg_def, Ideal.Quotient.lift_mk]
      show Ideal.Quotient.mk (Ideal.span (Set.range q)) (MvPolynomial.aeval ev' (MvPolynomial.C a)) = _
      rw [MvPolynomial.aeval_C]
      rfl
    · intro j
      show g (f (Ideal.Quotient.mk (Ideal.span (Set.range q)) (MvPolynomial.X j))) =
        Ideal.Quotient.mk (Ideal.span (Set.range q)) (MvPolynomial.X j)
      rw [hf_def, Ideal.Quotient.lift_mk]
      show g (Ideal.Quotient.mk (Ideal.span (Set.range q₂)) (MvPolynomial.aeval ev (MvPolynomial.X j))) = _
      rw [MvPolynomial.aeval_X, hg_def, Ideal.Quotient.lift_mk]
      show Ideal.Quotient.mk (Ideal.span (Set.range q)) (MvPolynomial.aeval ev' (ev j)) = _
      rw [Ideal.Quotient.eq]
      exact hround1 j

open scoped TensorProduct in
/-- **「項目(d)の第二段」の実用形——遷移写像の候補が実際に`RingEquiv`
として手に入る**。`exists_mvPolynomial_quotient_ringEquiv_descend`
(生データ)を`exists_mvPolynomial_quotient_ringEquiv_of_data`(生データ
から`RingEquiv`を組み立てる)へ渡すだけ——`Lemma 4.1`の`GlueData`構成
(`Spec`を取れば2つの`descendPieceR`片の間の実際のスキーム同型になる)
で直接使う最終形。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_mvPolynomial_quotient_ringEquiv_descend' (A : Type) [CommRing A] [Algebra ℚ A]
    (R R₂ : FgSubalgebra ℚ ℝ) {ι κ ι' κ' : Type} [Fintype ι] [Fintype ι'] [Fintype κ] [Fintype κ']
    (q : κ → MvPolynomial ι (A ⊗[ℚ] R.1)) (q₂ : κ' → MvPolynomial ι' (A ⊗[ℚ] R₂.1))
    (ψ : ι → MvPolynomial ι' (A ⊗[ℚ] ℝ)) (ψ' : ι' → MvPolynomial ι (A ⊗[ℚ] ℝ))
    (hψ : ∀ k, MvPolynomial.aeval ψ
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom (q k)) ∈
      Ideal.span (Set.range (fun k' => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R₂.1)).toRingHom (q₂ k'))))
    (hψ' : ∀ k', MvPolynomial.aeval ψ'
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R₂.1)).toRingHom (q₂ k')) ∈
      Ideal.span (Set.range (fun k => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom (q k))))
    (hround1 : ∀ i, MvPolynomial.aeval ψ' (ψ i) - MvPolynomial.X i ∈
      Ideal.span (Set.range (fun k => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom (q k))))
    (hround2 : ∀ j, MvPolynomial.aeval ψ (ψ' j) - MvPolynomial.X j ∈
      Ideal.span (Set.range (fun k' => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R₂.1)).toRingHom (q₂ k')))) :
    ∃ (R' : FgSubalgebra ℚ ℝ) (hR : R ≤ R') (hR₂ : R₂ ≤ R'),
      Nonempty ((MvPolynomial ι (A ⊗[ℚ] R'.1) ⧸ Ideal.span (Set.range (fun k => MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR)).toRingHom (q k)))) ≃+*
        (MvPolynomial ι' (A ⊗[ℚ] R'.1) ⧸ Ideal.span (Set.range (fun k' => MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR₂)).toRingHom (q₂ k'))))) := by
  obtain ⟨R', hR, hR₂, ev, ev', hev, hev', hq, hq₂, hround1', hround2'⟩ :=
    exists_mvPolynomial_quotient_ringEquiv_descend A R R₂ q q₂ ψ ψ' hψ hψ' hround1 hround2
  exact ⟨R', hR, hR₂, exists_mvPolynomial_quotient_ringEquiv_of_data _ _ ev ev' hq hq₂ hround1' hround2'⟩

open scoped TensorProduct in
open CategoryTheory AlgebraicGeometry in
/-- **「項目(d)の第二段」のスキーム版——遷移写像の候補が実際のスキーム
同型として手に入る**。`exists_mvPolynomial_quotient_ringEquiv_descend'`
(環レベルの`RingEquiv`)を`Scheme.Spec.mapIso`で送るだけ——`descendPieceR`
(`ExtLimit.lean`、`Spec(S_0)`の形)2つの間の実際のスキーム同型を、
共通の精密化`R'`上で構成的に得る、`Lemma 4.1`の`GlueData`構成で直接
使う最終形。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_mvPolynomial_quotient_specIso_descend (A : Type) [CommRing A] [Algebra ℚ A]
    (R R₂ : FgSubalgebra ℚ ℝ) {ι κ ι' κ' : Type} [Fintype ι] [Fintype ι'] [Fintype κ] [Fintype κ']
    (q : κ → MvPolynomial ι (A ⊗[ℚ] R.1)) (q₂ : κ' → MvPolynomial ι' (A ⊗[ℚ] R₂.1))
    (ψ : ι → MvPolynomial ι' (A ⊗[ℚ] ℝ)) (ψ' : ι' → MvPolynomial ι (A ⊗[ℚ] ℝ))
    (hψ : ∀ k, MvPolynomial.aeval ψ
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom (q k)) ∈
      Ideal.span (Set.range (fun k' => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R₂.1)).toRingHom (q₂ k'))))
    (hψ' : ∀ k', MvPolynomial.aeval ψ'
        (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R₂.1)).toRingHom (q₂ k')) ∈
      Ideal.span (Set.range (fun k => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom (q k))))
    (hround1 : ∀ i, MvPolynomial.aeval ψ' (ψ i) - MvPolynomial.X i ∈
      Ideal.span (Set.range (fun k => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom (q k))))
    (hround2 : ∀ j, MvPolynomial.aeval ψ (ψ' j) - MvPolynomial.X j ∈
      Ideal.span (Set.range (fun k' => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R₂.1)).toRingHom (q₂ k')))) :
    ∃ (R' : FgSubalgebra ℚ ℝ) (hR : R ≤ R') (hR₂ : R₂ ≤ R'),
      letI : CommRing (A ⊗[ℚ] R'.1) := inferInstance
      Nonempty (Spec (CommRingCat.of (MvPolynomial ι (A ⊗[ℚ] R'.1) ⧸ Ideal.span (Set.range (fun k =>
          MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR)).toRingHom
            (q k))))) ≅
        Spec (CommRingCat.of (MvPolynomial ι' (A ⊗[ℚ] R'.1) ⧸ Ideal.span (Set.range (fun k' =>
          MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR₂)).toRingHom
            (q₂ k')))))) := by
  obtain ⟨R', hR, hR₂, e⟩ := exists_mvPolynomial_quotient_ringEquiv_descend' A R R₂ q q₂ ψ ψ' hψ hψ' hround1 hround2
  letI : CommRing (A ⊗[ℚ] R'.1) := inferInstance
  exact ⟨R', hR, hR₂, ⟨Scheme.Spec.mapIso e.some.symm.toCommRingCatIso.op⟩⟩

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

/-! ## `descendPieceR`(`ExtLimit.lean`)を実際の局所化として貼り合わせる
ための代数的な核心事実(`2026-09-05続き17`)

`ExtLimit.lean`の`piece_basicOpen_mul_eq`で確立した「`C`側の片
`piece(D(f*g))`は`piece(D(f))`の基本開そのもの」という事実を、
`descendPieceR`(`R`レベルの抽象スキーム)側の**実際の開埋め込み**へ
橋渡しするための、`CorrHyp`非依存の一般的な可換環論の事実——
「`S₀`を`h`で局所化してから`B`上で`T`とテンソルする」のと「先に`B`上で
`T`とテンソルしてから、その中の`h`の像で局所化する」のが一致する
(`Localization.Away`は底変換と可換)。 -/

open scoped TensorProduct in
/-- **`Localization.Away`は底変換と可換**——`M`が`S₀`の`h`による局所化
(`IsLocalization.Away h M`)のとき、`M`を`B`上で`T`とテンソルしたもの
(`M ⊗[B] T`)は、先に`S₀`自体を`B`上で`T`とテンソルしてから、その中の
`h`の像であらためて局所化したもの(`Localization.Away (algebraMap S₀
(S₀⊗[B]T) h)`)に等しい。証明は`Algebra.TensorProduct.cancelBaseChange`
(結合律・簡約)+`IsLocalization.Away.tensorRight`(局所化とテンソルの
可換性、mathlibのインスタンス)+`IsLocalization.algEquiv`(同じ部分集合
による局所化の一意性)を合成するだけ——`Algebra.TensorProduct.
rightAlgebra`は`abbrev`(自動探索されない)なので`letI`で明示的に
供給する必要があった、新しい配管の一手。

`descendPieceR`を`W := X.basicOpen(f*g)`について「`D(f)`側の`descendPieceR`
の環をある元で局所化したもの」として**直接構成**すれば、この事実により
そのℝへの底変換が正しく`Γ(C,piece(W))`に一致することが従い、かつ
`Spec`への言い換え(`Spec.map (algebraMap S₀ M)`)は`IsOpenImmersion.
of_isLocalization`(mathlib)により**自動的に開埋め込みになる**——
`Lemma 4.1`のGlueDataで要求される`f i j : V(i,j) ⟶ U i`の開埋め込み性を、
独立に選んだ2つの`Presentation`を事後的に比較するのではなく**構成に
よって保証する**という設計転換の核心部品。

★**sorry 無し**。標準3公理のみ。 -/
theorem isLocalization_away_tensor_eq (B S₀ T M : Type) [CommRing B] [CommRing S₀] [CommRing T] [CommRing M]
    [Algebra B S₀] [Algebra B T] [Algebra B M] [Algebra S₀ M] [IsScalarTower B S₀ M]
    (h : S₀) [IsLocalization.Away h M] :
    Nonempty (M ⊗[B] T ≃+* Localization.Away (algebraMap S₀ (S₀ ⊗[B] T) h)) := by
  letI : Algebra (S₀ ⊗[B] T) (M ⊗[S₀] (S₀ ⊗[B] T)) := Algebra.TensorProduct.rightAlgebra
  have e1 : M ⊗[S₀] (S₀ ⊗[B] T) ≃ₐ[S₀] M ⊗[B] T := Algebra.TensorProduct.cancelBaseChange B S₀ S₀ M T
  haveI : IsLocalization.Away (algebraMap S₀ (S₀ ⊗[B] T) h) (M ⊗[S₀] (S₀ ⊗[B] T)) :=
    IsLocalization.Away.tensorRight (S := S₀ ⊗[B] T) h M
  have e2 := IsLocalization.algEquiv (Submonoid.powers (algebraMap S₀ (S₀ ⊗[B] T) h))
    (M ⊗[S₀] (S₀ ⊗[B] T)) (Localization.Away (algebraMap S₀ (S₀ ⊗[B] T) h))
  exact ⟨e1.toRingEquiv.symm.trans e2.toRingEquiv⟩

/-! ## `Γ(C,piece(D(f)))`の任意の元を`R`レベルへ持ち上げる(`2026-09-05続き18`)

`isLocalization_away_tensor_eq`を実際に使う(`descendPieceR`をD(f*g)`に
ついて`D(f)`側の局所化として直接構成する)には、局所化パラメータ
`piece_basicOpen_localizationElem`(`Γ(C,piece(D(f)))`の元、`ExtLimit.lean`)
自体を`descendPieceR`の`R`レベルの環`S₀`の元として持ち上げる必要がある
——`S₀⊗[B](A⊗ℝ) ≅ Γ(C,piece(D(f)))`(`pieceAlgebra_R_model_baseChange`)
の**逆像**を、`quotient_mvPolynomial_baseChange`(`S₀⊗[B](A⊗ℝ)`を
`MvPolynomial(A⊗ℝ)`商として実現)+`exists_fg_subalgebra_tensor_
mvPolynomial_finset`(既存、多項式の係数を`R`レベルへ降ろす)で構成する
——`pieceAlgebra_relation_descend_q₀`が「関係式の族」を降ろしたのと
同じ技法を、「単一の任意の元」に適用しただけ。 -/

open scoped TensorProduct Classical in
/-- **商`MvPolynomial`商⊗ℝの任意の元は、`R`を昇格すれば`MvPolynomial`
係数として`R`レベルへ持ち上げられる**——`Ideal.Quotient.mk_surjective`
で多項式の代表元を取り、その有限個の係数を`exists_fg_subalgebra_tensor_
mvPolynomial_finset`で降ろすだけ。`descendPieceR`の`S₀`(`MvPolynomial
(Fin n)(A⊗R.1)⧸I`型)の元を、`Γ(C,piece)`側の任意の元から実際に構成する
ための核心部品——`isLocalization_away_tensor_eq`の`h`をこれで得る。

配管の注意(新しい失敗形、`tools/lean-idioms.md`に追記の価値あり):
`MvPolynomial ι (A⊗[ℚ]R.1) ⧸ I`という「`TensorProduct`係数の
`MvPolynomial`の商」を書くと、`HasQuotient`の自動探索が(`Ideal`の
`Semiring`インスタンスを`AddMonoidAlgebra.semiring`経由で解決しようと
して)`MvPolynomial`本来の`CommRing`インスタンスと非`defeq`に見える形で
失敗する——`descendPieceR`と同じく`letI hCR : CommRing (A⊗[ℚ]R.1) :=
inferInstance`で**係数環自身**のインスタンスを先に確定させておくと解消
する(`MvPolynomial`側を先に確定させても効果が無い、係数環側が鍵)。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_fg_subalgebra_tensor_quotientMvPolynomial_lift (A : Type) [CommRing A] [Algebra ℚ A]
    (R : FgSubalgebra ℚ ℝ) {ι : Type} [Fintype ι] :
    letI hCR : CommRing (A ⊗[ℚ] R.1) := inferInstance
    letI : Algebra (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ) :=
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom.toAlgebra
    ∀ (I : Ideal (MvPolynomial ι (A ⊗[ℚ] R.1)))
      (z : (MvPolynomial ι (A ⊗[ℚ] R.1) ⧸ I) ⊗[A ⊗[ℚ] R.1] (A ⊗[ℚ] ℝ)),
    ∃ (R' : FgSubalgebra ℚ ℝ) (_hR : R ≤ R') (p₀ : MvPolynomial ι (A ⊗[ℚ] R'.1)),
      (quotient_mvPolynomial_baseChange (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ) ι I) z =
        Ideal.Quotient.mk (Ideal.map (MvPolynomial.map (algebraMap (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ))) I)
          (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom p₀) := by
  letI hCR : CommRing (A ⊗[ℚ] R.1) := inferInstance
  letI : Algebra (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ) :=
    (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom.toAlgebra
  intro I z
  obtain ⟨p, hp⟩ := Ideal.Quotient.mk_surjective
    ((quotient_mvPolynomial_baseChange (A ⊗[ℚ] R.1) (A ⊗[ℚ] ℝ) ι I) z)
  obtain ⟨R'', hR''⟩ := exists_fg_subalgebra_tensor_mvPolynomial_finset A ({p} : Finset (MvPolynomial ι (A ⊗[ℚ] ℝ)))
  obtain ⟨p₀, hp₀⟩ := hR'' p (Finset.mem_singleton_self p)
  obtain ⟨R', hRR', hR''R'⟩ := exists_fgSubalgebra_upperBound2 R R''
  refine ⟨R', hRR', MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A)
    (Subalgebra.inclusion hR''R')).toRingHom p₀, ?_⟩
  rw [← hp, ← hp₀, MvPolynomial.map_map]
  have h1 : (Subalgebra.val R'.1).comp (Subalgebra.inclusion hR''R') = Subalgebra.val R''.1 := rfl
  have h2 := Algebra.TensorProduct.map_comp (S := ℚ) (R := ℚ) (AlgHom.id ℚ A) (AlgHom.id ℚ A)
    (Subalgebra.val R'.1) (Subalgebra.inclusion hR''R')
  rw [AlgHom.id_comp, h1] at h2
  have h3 : (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R'.1)).toRingHom.comp
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.inclusion hR''R')).toRingHom =
      (Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R''.1)).toRingHom :=
    congrArg AlgHom.toRingHom h2.symm
  rw [h3]

/-- **`I`を`B→B'`で降ろしてから`B'→T`へ底変換したものは、`B→T`へ直接
底変換したものと一致する**——`Ideal.map_map`(合成)+`IsScalarTower.
algebraMap_eq`(`B→B'→T`が`B→T`と一致すること)を`MvPolynomial.
ringHom_ext`で確認するだけの、`CorrHyp`非依存の一般的な可換環論の事実。
`descendPieceR`(`ExtLimit.lean`)を`R`レベルからさらに`R'`レベルへ昇格
したものが、`descendPieceR`自体の`ℝ`への底変換(`pieceAlgebra_R_model_
baseChange`)と`quotient_mvPolynomial_baseChange`を通じて正しく整合する
ことを示すための橋渡し——`descendPieceR`の局所化(`Lemma 4.1`)がℝレベル
で正しい対象を実現していることの証明の核心部品。

★**sorry 無し**。標準3公理のみ。 -/
theorem ideal_map_mvPolynomial_promote_baseChange_eq (B B' T ι : Type) [CommRing B] [CommRing B'] [CommRing T]
    [Algebra B B'] [Algebra B' T] [Algebra B T] [IsScalarTower B B' T] (I : Ideal (MvPolynomial ι B)) :
    Ideal.map (MvPolynomial.map (algebraMap B' T)) (Ideal.map (MvPolynomial.map (algebraMap B B')) I) =
      Ideal.map (MvPolynomial.map (algebraMap B T)) I := by
  rw [Ideal.map_map]
  congr 1
  apply MvPolynomial.ringHom_ext
  · intro b
    simp [RingHom.comp_apply, IsScalarTower.algebraMap_eq B B' T]
  · intro x
    simp

/-- **環同型が局所化元同士を対応させれば、対応する`Away`局所化も
同型になる**——`IsLocalization.ringEquivOfRingEquiv`(mathlib)を
直接適用するだけだが、`A`・`B`を**抽象的な`CommRing`**として一般化
した形で用意しておくことが鍵——`A`・`B`をテンソル積などの具体的な
構成のまま直接この補題を書こうとすると、`Mul`/`Add`インスタンスが
(`Algebra.TensorProduct.instMul`経由か`instDistribOfSemiring.toMul`
経由かで)非`defeq`に見えるinstance diamondに当たる(`descendPieceR`の
局所化を実際にℝレベルへ橋渡しする際に発見した、新しい配管の教訓)——
この補題のように**一度抽象的な`CommRing`型として切り出してから**
具体的な型を代入すれば、そのdiamondを回避できる。

★**sorry 無し**。標準3公理のみ。 -/
theorem ringEquiv_localization_of_apply_eq (A B : Type) [CommRing A] [CommRing B] (e : A ≃+* B) (a : A) (b : B)
    (hab : e a = b) : Nonempty (Localization.Away a ≃+* Localization.Away b) :=
  ⟨IsLocalization.ringEquivOfRingEquiv (Localization.Away a) (Localization.Away b) e
    (hab ▸ Submonoid.map_powers e.toMonoidHom a)⟩

open scoped TensorProduct in
/-- **`descendPieceR`の`R'`レベル局所化は、実際にℝレベルで正しい対象
(`quotient_mvPolynomial_baseChange`の右辺の局所化)を実現する**——
`Lemma 4.1`の`GlueData`構築で残っていた最後の橋渡し。`isLocalization_
away_tensor_eq`(局所化と底変換の可換性)・`ideal_map_mvPolynomial_
promote_baseChange_eq`(イデアルの2段底変換の一致)・`quotient_
mvPolynomial_baseChange_tmul_one`(純テンソル上の自然性)・
`ringEquiv_localization_of_apply_eq`(局所化元の対応から局所化自身の
同型を得る)という4つの部品を繋いで完成させた。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_ringEquiv_localization_of_eq (B B' T ι : Type) [CommRing B] [CommRing B'] [CommRing T]
    [Algebra B B'] [Algebra B' T] [Algebra B T] [IsScalarTower B B' T]
    (I : Ideal (MvPolynomial ι B)) (p₀ : MvPolynomial ι B') :
    letI hIR' : Ideal (MvPolynomial ι B') := Ideal.map (MvPolynomial.map (algebraMap B B')) I
    ∀ (M : Type) [CommRing M] [Algebra B' M]
      [Algebra (MvPolynomial ι B' ⧸ hIR') M]
      [IsScalarTower B' (MvPolynomial ι B' ⧸ hIR') M]
      [IsLocalization.Away (Ideal.Quotient.mk hIR' p₀) M],
    Nonempty (M ⊗[B'] T ≃+* Localization.Away (Ideal.Quotient.mk
      (Ideal.map (MvPolynomial.map (algebraMap B T)) I) (MvPolynomial.map (algebraMap B' T) p₀))) := by
  intro M _ _ _ _ _
  obtain ⟨e0⟩ := isLocalization_away_tensor_eq B' (MvPolynomial ι B' ⧸
      Ideal.map (MvPolynomial.map (algebraMap B B')) I) T M
    (Ideal.Quotient.mk (Ideal.map (MvPolynomial.map (algebraMap B B')) I) p₀)
  rw [Algebra.TensorProduct.algebraMap_apply] at e0
  set IR' := Ideal.map (MvPolynomial.map (algebraMap B B')) I with hIR'def
  have hIdealEq : Ideal.map (MvPolynomial.map (algebraMap B' T)) IR' =
      Ideal.map (MvPolynomial.map (algebraMap B T)) I :=
    ideal_map_mvPolynomial_promote_baseChange_eq B B' T ι I
  set e2 : (MvPolynomial ι B' ⧸ IR') ⊗[B'] T ≃+* MvPolynomial ι T ⧸ Ideal.map (MvPolynomial.map (algebraMap B T)) I :=
    (quotient_mvPolynomial_baseChange B' T ι IR').trans (Ideal.quotEquivOfEq hIdealEq) with he2def
  have key : e2 ((Ideal.Quotient.mk IR' p₀) ⊗ₜ[B'] (1 : T)) =
      Ideal.Quotient.mk (Ideal.map (MvPolynomial.map (algebraMap B T)) I) (MvPolynomial.map (algebraMap B' T) p₀) := by
    rw [he2def, RingEquiv.trans_apply, quotient_mvPolynomial_baseChange_tmul_one, Ideal.quotEquivOfEq_mk]
  obtain ⟨e3⟩ := ringEquiv_localization_of_apply_eq ((MvPolynomial ι B' ⧸ IR') ⊗[B'] T)
    (MvPolynomial ι T ⧸ Ideal.map (MvPolynomial.map (algebraMap B T)) I) e2
    ((Ideal.Quotient.mk IR' p₀) ⊗ₜ[B'] (1 : T))
    (Ideal.Quotient.mk (Ideal.map (MvPolynomial.map (algebraMap B T)) I) (MvPolynomial.map (algebraMap B' T) p₀))
    key
  exact ⟨e0.trans e3⟩

open scoped TensorProduct Classical in
/-- **`descendPieceR`の`R'`レベル局所化は、`Γ(C,piece)`のアフィン片を
局所化して得られる任意の対象(`Wtarget`、`piece(D(f*g))`の実現)とも
同型になる**——`exists_ringEquiv_localization_of_eq`(`R`↔ℝ橋渡し)・
`ringEquiv_localization_of_apply_eq`(局所化元の対応)・`IsLocalization.
algEquiv`(局所化の実現の一意性)の3部品を、`hp₀`(`exists_piece_
basicOpen_R_lift`が与える対応)を仲立ちに繋ぐだけ。`CorrHyp`固有の
`pieceAlgebra`等の足場を一切使わない**抽象的な**形にしたことで、
巨大な足場を伴う証明で繰り返し当たっていたelaboration timeoutを完全に
回避できた——呼び出し側は`e`・`h₂`・`hp₀`・`Wtarget`の`IsLocalization.
Away`インスタンスを揃えるだけで、この補題を直接適用できる。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_ringEquiv_of_piece_lift (B B' T : Type) [CommRing B] [CommRing B'] [CommRing T]
    [Algebra B B'] [Algebra B' T] [Algebra B T] [IsScalarTower B B' T]
    (n : Type) [Fintype n] (I : Ideal (MvPolynomial n B)) (p₀ : MvPolynomial n B')
    (Ctarget Wtarget : Type) [CommRing Ctarget] [CommRing Wtarget]
    (e : (MvPolynomial n B ⧸ I) ⊗[B] T ≃+* Ctarget) (h₂ : Ctarget)
    [Algebra Ctarget Wtarget] [IsLocalization.Away h₂ Wtarget]
    (hp₀ : (quotient_mvPolynomial_baseChange B T n I) (e.symm h₂) =
      Ideal.Quotient.mk (Ideal.map (MvPolynomial.map (algebraMap B T)) I) (MvPolynomial.map (algebraMap B' T) p₀)) :
    letI hIR' : Ideal (MvPolynomial n B') := Ideal.map (MvPolynomial.map (algebraMap B B')) I
    ∀ (M : Type) [CommRing M] [Algebra B' M]
      [Algebra (MvPolynomial n B' ⧸ hIR') M]
      [IsScalarTower B' (MvPolynomial n B' ⧸ hIR') M]
      [IsLocalization.Away (Ideal.Quotient.mk hIR' p₀) M],
    Nonempty (M ⊗[B'] T ≃+* Wtarget) := by
  intro M _ _ _ _ _
  obtain ⟨e4⟩ := exists_ringEquiv_localization_of_eq B B' T n I p₀ M
  obtain ⟨e5⟩ := ringEquiv_localization_of_apply_eq _ _ e.symm h₂ (e.symm h₂) rfl
  obtain ⟨e6⟩ := ringEquiv_localization_of_apply_eq _ _
    (quotient_mvPolynomial_baseChange B T n I) (e.symm h₂) _ hp₀
  have e7 := IsLocalization.algEquiv (Submonoid.powers h₂) Wtarget (Localization.Away h₂)
  exact ⟨e4.trans (e6.symm.trans (e5.symm.trans e7.toRingEquiv.symm))⟩

/-! ## `t`(`Lemma 4.1`の`GlueData`の遷移射)へ向けた第一歩——
`Localization.Away`を`MvPolynomial`の商として実現する
(`2026-09-05夜、続き10`)

`M_ij`(`D(f)`側)・`M_ji`(`D(g)`側)を`exists_mvPolynomial_quotient_
specIso_descend`(既存、`ψ・ψ'`による2つの独立な有限表示の比較)へ
渡すには、両方を**`MvPolynomial`の関係式の族**として具体的に表す
必要がある——`Localization.Away h`(`h`は`MvPolynomial n B⧸I`の元)の
標準的な「逆元を1変数添加する」表示を確立する。 -/

open scoped Classical in
/-- **`MvPolynomial n B⧸I`の元`p`(のクラス)による`Away`局所化は、
`MvPolynomial Unit(MvPolynomial n B)`を`I`の像と局所化関係式で割った
商そのものと同一視できる**——`IsLocalization.Away.mvPolynomialQuotient
Equiv`(mathlib、`S:=Localization.Away h`の場合に特殊化)と`MvPolynomial.
quotientEquivQuotientMvPolynomial`(mathlib、`MvPolynomial σ(R⧸I)≃
MvPolynomial σ R⧸(Iの像)`)を、両者が局所化関係式の定義元(`C(mk I p)`・
`X()`)を同じ場所へ送ること(`hnat`・`hnat2`)を確認しながら`Ideal.
quotientEquiv`で繋ぐ。`t`の構成(2つの独立な`R`レベル局所化を`MvPolynomial`
の関係式の族として比較する)への第一歩——次の一手は、この「入れ子の
商」を`MvPolynomial.sumAlgEquiv`で`MvPolynomial(n⊕Unit)B`という**1段の
商**へ平坦化すること(`exists_mvPolynomial_quotient_specIso_descend`の
`q`・`q₂`が要求する形)。

★**sorry 無し**。標準3公理のみ。 -/
theorem localization_away_quotient_mvPolynomial_equiv (B : Type) [CommRing B] (n : Type) [Fintype n]
    (I : Ideal (MvPolynomial n B)) (p : MvPolynomial n B) :
    Nonempty (Localization.Away (Ideal.Quotient.mk I p) ≃+*
      (MvPolynomial Unit (MvPolynomial n B) ⧸ Ideal.map MvPolynomial.C I) ⧸
        (Ideal.span {MvPolynomial.C p * MvPolynomial.X () - 1}).map (Ideal.Quotient.mk (Ideal.map MvPolynomial.C I))) := by
  set e1 := IsLocalization.Away.mvPolynomialQuotientEquiv (Localization.Away (Ideal.Quotient.mk I p))
    (Ideal.Quotient.mk I p) with he1
  set e2 := MvPolynomial.quotientEquivQuotientMvPolynomial (σ := Unit) I with he2
  have hnat : e2 (MvPolynomial.C (Ideal.Quotient.mk I p)) =
      Ideal.Quotient.mk (Ideal.map MvPolynomial.C I) (MvPolynomial.C p) := by
    simp [he2, MvPolynomial.quotientEquivQuotientMvPolynomial]
  have hnat2 : e2 (MvPolynomial.X ()) = Ideal.Quotient.mk (Ideal.map MvPolynomial.C I) (MvPolynomial.X ()) := by
    simp [he2, MvPolynomial.quotientEquivQuotientMvPolynomial]
  have hIdealEq : (Ideal.span {MvPolynomial.C p * MvPolynomial.X () - 1}).map (Ideal.Quotient.mk (Ideal.map MvPolynomial.C I)) =
      (Ideal.span {MvPolynomial.C (Ideal.Quotient.mk I p) * MvPolynomial.X () - 1}).map e2.toRingEquiv.toRingHom := by
    rw [Ideal.map_span, Ideal.map_span, Set.image_singleton, Set.image_singleton]
    have hval1 : Ideal.Quotient.mk (Ideal.map MvPolynomial.C I) (MvPolynomial.C p * MvPolynomial.X () - 1) =
        Ideal.Quotient.mk (Ideal.map MvPolynomial.C I) (MvPolynomial.C p) *
          Ideal.Quotient.mk (Ideal.map MvPolynomial.C I) (MvPolynomial.X ()) - 1 := by
      rw [map_sub, map_mul, map_one]
    have hval2 : e2.toRingEquiv.toRingHom (MvPolynomial.C (Ideal.Quotient.mk I p) * MvPolynomial.X () - 1) =
        Ideal.Quotient.mk (Ideal.map MvPolynomial.C I) (MvPolynomial.C p) *
          Ideal.Quotient.mk (Ideal.map MvPolynomial.C I) (MvPolynomial.X ()) - 1 := by
      show e2 (MvPolynomial.C (Ideal.Quotient.mk I p) * MvPolynomial.X () - 1) = _
      rw [map_sub, map_mul, map_one, hnat, hnat2]
    rw [hval1, hval2]
  exact ⟨e1.symm.toRingEquiv.trans (Ideal.quotientEquiv _ _ e2.toRingEquiv hIdealEq)⟩

open scoped Classical in
/-- `localization_away_quotient_mvPolynomial_equiv`の入れ子の商
(`MvPolynomial Unit(MvPolynomial n B)`の商をさらに商)を、
`exists_mvPolynomial_quotient_specIso_descend`が要求する**1段の**
`MvPolynomial(Unit⊕n)B`の商へ平坦化する。`DoubleQuot.quotQuotEquivQuotSup`
(mathlib)で入れ子の商を`I'⊔span{...}`単独の商へ、`MvPolynomial.sumAlgEquiv`
(mathlib)で索引を`Unit⊕n`単独へ、それぞれ変換する。生成元の対応は
`sumAlgEquiv_comp_rename_inl`・`_inr`(いずれもmathlib、naturality)から
具体的に計算できる:
`e2.symm(C x) = rename Sum.inr x`、`e2.symm(X()) = X(Sum.inl())`。

★**sorry 無し**。標準3公理のみ。 -/
theorem localization_away_quotient_mvPolynomial_flat_equiv (B : Type) [CommRing B]
    (n : Type) [Fintype n] (κ₀ : Type) [Fintype κ₀]
    (q₀ : κ₀ → MvPolynomial n B) (p : MvPolynomial n B) :
    Nonempty (Localization.Away (Ideal.Quotient.mk (Ideal.span (Set.range q₀)) p) ≃+*
      MvPolynomial (Unit ⊕ n) B ⧸ Ideal.span (Set.range
        (Sum.elim (fun k => MvPolynomial.rename Sum.inr (q₀ k))
          (fun _ : Unit => MvPolynomial.rename Sum.inr p * MvPolynomial.X (Sum.inl ()) - 1)))) := by
  set I := Ideal.span (Set.range q₀) with hI
  obtain ⟨e0⟩ := localization_away_quotient_mvPolynomial_equiv B n I p
  set I' := Ideal.map MvPolynomial.C I with hI'def
  set e1 := DoubleQuot.quotQuotEquivQuotSup I' (Ideal.span {MvPolynomial.C p * MvPolynomial.X () - 1}) with he1
  set e2 := MvPolynomial.sumAlgEquiv B Unit n with he2
  have hI'span : I' = Ideal.span (Set.range (fun k => MvPolynomial.C (q₀ k))) := by
    rw [hI'def, hI, Ideal.map_span]
    congr 1
    exact (Set.range_comp _ _).symm
  have hJ : I' ⊔ Ideal.span {MvPolynomial.C p * MvPolynomial.X () - 1} =
      Ideal.span (Set.range (fun k => MvPolynomial.C (q₀ k)) ∪ {MvPolynomial.C p * MvPolynomial.X () - 1}) := by
    rw [hI'span, Ideal.span_union]
  have he2symC : ∀ x : MvPolynomial n B, e2.symm (MvPolynomial.C x) = MvPolynomial.rename Sum.inr x := by
    intro x
    have h := congrFun (congrArg DFunLike.coe (MvPolynomial.sumAlgEquiv_comp_rename_inr B Unit n)) x
    simp only [AlgHom.comp_apply, AlgEquiv.coe_algHom, IsScalarTower.coe_toAlgHom',
      MvPolynomial.algebraMap_eq] at h
    rw [← h, he2, AlgEquiv.symm_apply_apply]
  have he2symX : e2.symm (MvPolynomial.X () : MvPolynomial Unit (MvPolynomial n B)) = MvPolynomial.X (Sum.inl ()) := by
    have h := congrFun (congrArg DFunLike.coe (MvPolynomial.sumAlgEquiv_comp_rename_inl B Unit n)) (MvPolynomial.X ())
    simp only [AlgHom.comp_apply, AlgEquiv.coe_algHom, MvPolynomial.mapAlgHom_apply, MvPolynomial.map_X,
      MvPolynomial.rename_X] at h
    rw [← h, he2, AlgEquiv.symm_apply_apply]
  have hEq : Ideal.span (Set.range
        (Sum.elim (fun k => MvPolynomial.rename Sum.inr (q₀ k))
          (fun _ : Unit => MvPolynomial.rename Sum.inr p * MvPolynomial.X (Sum.inl ()) - 1))) =
      Ideal.map e2.symm.toRingEquiv.toRingHom
        (I' ⊔ Ideal.span {MvPolynomial.C p * MvPolynomial.X () - 1}) := by
    rw [hJ, Ideal.map_span, Set.Sum.elim_range]
    congr 1
    rw [Set.image_union, Set.image_singleton, Set.range_const]
    have himg1 : (e2.symm.toRingEquiv.toRingHom : MvPolynomial Unit (MvPolynomial n B) →+* MvPolynomial (Unit ⊕ n) B) ''
        Set.range (fun k => MvPolynomial.C (q₀ k)) =
        Set.range (fun k => MvPolynomial.rename Sum.inr (q₀ k)) := by
      rw [← Set.range_comp]
      congr 1
      funext k
      exact he2symC (q₀ k)
    have himg2 : (e2.symm.toRingEquiv.toRingHom : MvPolynomial Unit (MvPolynomial n B) →+* MvPolynomial (Unit ⊕ n) B)
        (MvPolynomial.C p * MvPolynomial.X () - 1) =
        MvPolynomial.rename Sum.inr p * MvPolynomial.X (Sum.inl ()) - 1 := by
      show e2.symm (MvPolynomial.C p * MvPolynomial.X () - 1) = _
      rw [map_sub, map_mul, map_one, he2symC, he2symX]
    rw [himg1, himg2]
  exact ⟨e0.trans (e1.trans (Ideal.quotientEquiv _ _ e2.symm.toRingEquiv hEq))⟩

open scoped Classical in
/-- **`localization_away_quotient_mvPolynomial_flat_equiv`を、生成元の
係数が`φ : B→+*B'`で昇格された(`Ideal.map (MvPolynomial.map φ)`)場合へ
一般化**——`exists_descendPieceR_localization_baseChange`が実際に扱う
「`R`レベルの関係式`I`を`R'`レベルへ`algebraMap`で昇格したイデアル」
という形にそのまま適用できるようにするための橋渡し。`Ideal.map_span`
で`Ideal.map(map φ)(span(range q₀)) = span(range(map φ∘q₀))`と書き
換えるだけで`flat_equiv`本体に帰着する。

★**sorry 無し**。標準3公理のみ。 -/
theorem localization_away_quotient_mvPolynomial_flat_equiv_of_map
    (B B' : Type) [CommRing B] [CommRing B'] (φ : B →+* B')
    (κ₀ ι : Type) [Fintype κ₀] [Fintype ι] (q₀ : κ₀ → MvPolynomial ι B) (p : MvPolynomial ι B') :
    Nonempty (Localization.Away (Ideal.Quotient.mk
        (Ideal.map (MvPolynomial.map φ) (Ideal.span (Set.range q₀))) p) ≃+*
      MvPolynomial (Unit ⊕ ι) B' ⧸ Ideal.span (Set.range
        (Sum.elim (fun k => MvPolynomial.rename Sum.inr (MvPolynomial.map φ (q₀ k)))
          (fun _ : Unit => MvPolynomial.rename Sum.inr p * MvPolynomial.X (Sum.inl ()) - 1)))) := by
  have hmap : Ideal.map (MvPolynomial.map φ) (Ideal.span (Set.range q₀)) =
      Ideal.span (Set.range (fun k => MvPolynomial.map φ (q₀ k))) := by
    rw [Ideal.map_span]
    congr 1
    exact (Set.range_comp _ _).symm
  rw [hmap]
  exact localization_away_quotient_mvPolynomial_flat_equiv B' ι κ₀ (fun k => MvPolynomial.map φ (q₀ k)) p

/-- **`IsLocalization`のインスタンスは環同型に沿って移送できる**——`S`が
`R`の`M`による局所化なら、任意の環同型`e:S≃+*P`に対し、`e`を通じて
定義した`P`上の`R`代数構造(`algebraMap R P := e∘algebraMap R S`、
定義から自動的に`e`と両立する)のもとで`P`も同じ局所化になる。
`AlgEquiv.ofRingEquiv`(mathlib、両立条件から`RingEquiv`を`AlgEquiv`へ
格上げ)+`IsLocalization.isLocalization_of_algEquiv`(mathlib)を組み
合わせるだけ——`localization_away_quotient_mvPolynomial_flat_equiv_of_map`
が与える具体的な環同型を、`exists_descendPieceR_localization_baseChange`
が要求する`IsLocalization.Away`インスタンスへ変換するための橋渡し。

★**sorry 無し**。標準3公理のみ。 -/
theorem isLocalization_of_ringEquiv_transport (R S P : Type) [CommRing R] [CommRing S] [CommRing P]
    (M : Submonoid R) [Algebra R S] [IsLocalization M S] (e : S ≃+* P) :
    letI : Algebra R P := (e.toRingHom.comp (algebraMap R S)).toAlgebra
    IsLocalization M P := by
  letI : Algebra R P := (e.toRingHom.comp (algebraMap R S)).toAlgebra
  have he : ∀ r, e (algebraMap R S r) = algebraMap R P r := fun r => rfl
  exact IsLocalization.isLocalization_of_algEquiv M (AlgEquiv.ofRingEquiv (f := e) he)

/-! ## 平坦化の`AlgEquiv`版——`≃+*`ではなく`≃ₐ[B]`で作り直す
(`2026-09-05夜、続き14`)

`localization_away_quotient_mvPolynomial_equiv`系は`≃+*`(環同型)として
作ったが、これを`ExtLimit.lean`側で使おうとすると破綻する:商環
`MvPolynomial ι B ⧸ J`は係数環`B`に対する**自前の**`SMul`
(`Submodule.Quotient.instSMul'`)を持っており、環同型越しに移送した
`Algebra B _`はそれに負ける(`IsScalarTower`の型が
`Submodule.Quotient.instSMul'`で表示され、`of_algebraMap_eq`が返す
`Algebra.toSMul`3本組と合わない)。`tools/lean-idioms.md`の`#52`。

そこで**最初から係数環`B`上の`AlgEquiv`として作る**。mathlibの部品は
すべて`AlgEquiv`版が揃っている——`IsLocalization.Away.
mvPolynomialQuotientEquiv`(`≃ₐ[R]`)・`MvPolynomial.
quotientEquivQuotientMvPolynomial`(`≃ₐ[R]`)・`DoubleQuot.
quotQuotEquivQuotSupₐ`(`≃ₐ[R]`)・`Ideal.quotientEquivAlg`・
`Ideal.quotientEquivAlgOfEq`・`MvPolynomial.sumAlgEquiv`。基底が
途中で変わる箇所は`AlgEquiv.restrictScalars`で下の基底へ落とす。 -/

open scoped Classical in
/-- **`localization_away_quotient_mvPolynomial_equiv`の`AlgEquiv`版**——
`MvPolynomial n B⧸I`の元`p`による`Away`局所化の任意の実現`S`(係数環`B`
上のスカラー塔付き)が、`MvPolynomial Unit(MvPolynomial n B)`を`I`の像と
局所化関係式で割った商と`B`-代数として同型であること。証明は`≃+*`版と
同じ道筋(`mvPolynomialQuotientEquiv`+`quotientEquivQuotientMvPolynomial`
を`Ideal.quotientEquivAlg`で繋ぐ)を`restrictScalars B`で`B`上へ落として
行うだけ。

★**sorry 無し**。標準3公理のみ。 -/
theorem localization_away_quotient_mvPolynomial_algEquiv (B : Type) [CommRing B] (n : Type) [Fintype n]
    (I : Ideal (MvPolynomial n B)) (p : MvPolynomial n B)
    (S : Type) [CommRing S] [Algebra (MvPolynomial n B ⧸ I) S]
    [IsLocalization.Away (Ideal.Quotient.mk I p) S]
    [Algebra B S] [IsScalarTower B (MvPolynomial n B ⧸ I) S] :
    Nonempty (S ≃ₐ[B]
      (MvPolynomial Unit (MvPolynomial n B) ⧸ Ideal.map MvPolynomial.C I) ⧸
        (Ideal.span {MvPolynomial.C p * MvPolynomial.X () - 1}).map
          (Ideal.Quotient.mk (Ideal.map MvPolynomial.C I))) := by
  letI hCRQ : CommRing (MvPolynomial n B ⧸ I) := inferInstance
  haveI h2 : IsScalarTower B (MvPolynomial n B ⧸ I) (MvPolynomial Unit (MvPolynomial n B ⧸ I) ⧸
      Ideal.span {MvPolynomial.C (Ideal.Quotient.mk I p) * MvPolynomial.X () - 1}) := inferInstance
  set e1 : (MvPolynomial Unit (MvPolynomial n B ⧸ I) ⧸
      Ideal.span {MvPolynomial.C (Ideal.Quotient.mk I p) * MvPolynomial.X () - 1}) ≃ₐ[B] S :=
    (IsLocalization.Away.mvPolynomialQuotientEquiv S (Ideal.Quotient.mk I p)).restrictScalars B with he1
  set e2 : MvPolynomial Unit (MvPolynomial n B ⧸ I) ≃ₐ[B]
      (MvPolynomial Unit (MvPolynomial n B) ⧸ Ideal.map MvPolynomial.C I) :=
    (MvPolynomial.quotientEquivQuotientMvPolynomial (σ := Unit) I).restrictScalars B with he2
  have hnat : e2 (MvPolynomial.C (Ideal.Quotient.mk I p)) =
      Ideal.Quotient.mk (Ideal.map MvPolynomial.C I) (MvPolynomial.C p) := by
    simp [he2, MvPolynomial.quotientEquivQuotientMvPolynomial]
  have hnat2 : e2 (MvPolynomial.X ()) =
      Ideal.Quotient.mk (Ideal.map MvPolynomial.C I) (MvPolynomial.X ()) := by
    simp [he2, MvPolynomial.quotientEquivQuotientMvPolynomial]
  have hIdealEq : (Ideal.span {MvPolynomial.C p * MvPolynomial.X () - 1}).map
        (Ideal.Quotient.mk (Ideal.map MvPolynomial.C I)) =
      (Ideal.span {MvPolynomial.C (Ideal.Quotient.mk I p) * MvPolynomial.X () - 1}).map
        (e2 : MvPolynomial Unit (MvPolynomial n B ⧸ I) →ₐ[B]
          (MvPolynomial Unit (MvPolynomial n B) ⧸ Ideal.map MvPolynomial.C I)) := by
    rw [Ideal.map_span, Ideal.map_span, Set.image_singleton, Set.image_singleton]
    have hval1 : Ideal.Quotient.mk (Ideal.map MvPolynomial.C I)
        (MvPolynomial.C p * MvPolynomial.X () - 1) =
        Ideal.Quotient.mk (Ideal.map MvPolynomial.C I) (MvPolynomial.C p) *
          Ideal.Quotient.mk (Ideal.map MvPolynomial.C I) (MvPolynomial.X ()) - 1 := by
      rw [map_sub, map_mul, map_one]
    have hval2 : (e2 : MvPolynomial Unit (MvPolynomial n B ⧸ I) →ₐ[B] _)
        (MvPolynomial.C (Ideal.Quotient.mk I p) * MvPolynomial.X () - 1) =
        Ideal.Quotient.mk (Ideal.map MvPolynomial.C I) (MvPolynomial.C p) *
          Ideal.Quotient.mk (Ideal.map MvPolynomial.C I) (MvPolynomial.X ()) - 1 := by
      show e2 (MvPolynomial.C (Ideal.Quotient.mk I p) * MvPolynomial.X () - 1) = _
      rw [map_sub, map_mul, map_one, hnat, hnat2]
    rw [hval1, hval2]
  exact ⟨e1.symm.trans (Ideal.quotientEquivAlg _ _ e2 hIdealEq)⟩

open scoped Classical in
/-- **平坦化した1段の`MvPolynomial(Unit⊕n)B`商との`B`-代数同型**——
`localization_away_quotient_mvPolynomial_flat_equiv`の`AlgEquiv`版。
`DoubleQuot.quotQuotEquivQuotSupₐ`(`≃ₐ[B]`版、ただしイデアルが
`Ideal.map (Ideal.Quotient.mkₐ B I') J`という`mkₐ`の形なので
`Ideal.quotientEquivAlgOfEq`で`mk`の形へ橋渡しする)と
`MvPolynomial.sumAlgEquiv`(元から`≃ₐ[B]`)を使う。

配管の注意: `Ideal.map (Ideal.Quotient.mkₐ B I') J`の`IsTwoSided`
インスタンスは、係数環が入れ子の`MvPolynomial`のとき自動では見つからない
——`letI : CommRing (MvPolynomial n B)`と
`letI : CommRing (MvPolynomial Unit (MvPolynomial n B))`を先に登録する
(`tools/lean-idioms.md`の`#40`と同型)。

★**sorry 無し**。標準3公理のみ。 -/
theorem localization_away_quotient_mvPolynomial_flat_algEquiv (B : Type) [CommRing B]
    (n : Type) [Fintype n] (κ₀ : Type) [Fintype κ₀]
    (q₀ : κ₀ → MvPolynomial n B) (p : MvPolynomial n B)
    (S : Type) [CommRing S] [Algebra (MvPolynomial n B ⧸ Ideal.span (Set.range q₀)) S]
    [IsLocalization.Away (Ideal.Quotient.mk (Ideal.span (Set.range q₀)) p) S]
    [Algebra B S] [IsScalarTower B (MvPolynomial n B ⧸ Ideal.span (Set.range q₀)) S] :
    Nonempty (S ≃ₐ[B] MvPolynomial (Unit ⊕ n) B ⧸ Ideal.span (Set.range
        (Sum.elim (fun k => MvPolynomial.rename Sum.inr (q₀ k))
          (fun _ : Unit => MvPolynomial.rename Sum.inr p * MvPolynomial.X (Sum.inl ()) - 1)))) := by
  obtain ⟨e0⟩ := localization_away_quotient_mvPolynomial_algEquiv B n (Ideal.span (Set.range q₀)) p S
  letI hCRn : CommRing (MvPolynomial n B) := inferInstance
  letI hCRUn : CommRing (MvPolynomial Unit (MvPolynomial n B)) := inferInstance
  set I := Ideal.span (Set.range q₀) with hI
  set I' := Ideal.map MvPolynomial.C I with hI'def
  set J := Ideal.span {MvPolynomial.C p * MvPolynomial.X () - 1} with hJdef
  have hmkₐ : Ideal.map (Ideal.Quotient.mk I') J = Ideal.map (Ideal.Quotient.mkₐ B I') J := by
    rw [Ideal.map, Ideal.map, Ideal.Quotient.mkₐ_eq_mk]
  set e1 : ((MvPolynomial Unit (MvPolynomial n B) ⧸ I') ⧸ Ideal.map (Ideal.Quotient.mk I') J) ≃ₐ[B]
      MvPolynomial Unit (MvPolynomial n B) ⧸ I' ⊔ J :=
    (Ideal.quotientEquivAlgOfEq B hmkₐ).trans (DoubleQuot.quotQuotEquivQuotSupₐ B I' J) with he1
  set e2 := MvPolynomial.sumAlgEquiv B Unit n with he2
  have hI'span : I' = Ideal.span (Set.range (fun k => MvPolynomial.C (q₀ k))) := by
    rw [hI'def, hI, Ideal.map_span]
    congr 1
    exact (Set.range_comp _ _).symm
  have hJ : I' ⊔ J = Ideal.span (Set.range (fun k => MvPolynomial.C (q₀ k)) ∪
      {MvPolynomial.C p * MvPolynomial.X () - 1}) := by
    rw [hI'span, hJdef, Ideal.span_union]
  have he2symC : ∀ x : MvPolynomial n B, e2.symm (MvPolynomial.C x) = MvPolynomial.rename Sum.inr x := by
    intro x
    have h := congrFun (congrArg DFunLike.coe (MvPolynomial.sumAlgEquiv_comp_rename_inr B Unit n)) x
    simp only [AlgHom.comp_apply, AlgEquiv.coe_algHom, IsScalarTower.coe_toAlgHom',
      MvPolynomial.algebraMap_eq] at h
    rw [← h, he2, AlgEquiv.symm_apply_apply]
  have he2symX : e2.symm (MvPolynomial.X () : MvPolynomial Unit (MvPolynomial n B)) =
      MvPolynomial.X (Sum.inl ()) := by
    have h := congrFun (congrArg DFunLike.coe (MvPolynomial.sumAlgEquiv_comp_rename_inl B Unit n))
      (MvPolynomial.X ())
    simp only [AlgHom.comp_apply, AlgEquiv.coe_algHom, MvPolynomial.mapAlgHom_apply,
      MvPolynomial.map_X, MvPolynomial.rename_X] at h
    rw [← h, he2, AlgEquiv.symm_apply_apply]
  have hEq : Ideal.span (Set.range
        (Sum.elim (fun k => MvPolynomial.rename Sum.inr (q₀ k))
          (fun _ : Unit => MvPolynomial.rename Sum.inr p * MvPolynomial.X (Sum.inl ()) - 1))) =
      Ideal.map (e2.symm : MvPolynomial Unit (MvPolynomial n B) →ₐ[B] MvPolynomial (Unit ⊕ n) B)
        (I' ⊔ J) := by
    rw [hJ, Ideal.map_span, Set.Sum.elim_range]
    congr 1
    rw [Set.image_union, Set.image_singleton, Set.range_const]
    have himg1 : (e2.symm : MvPolynomial Unit (MvPolynomial n B) →ₐ[B] MvPolynomial (Unit ⊕ n) B) ''
        Set.range (fun k => MvPolynomial.C (q₀ k)) =
        Set.range (fun k => MvPolynomial.rename Sum.inr (q₀ k)) := by
      rw [← Set.range_comp]
      congr 1
      funext k
      exact he2symC (q₀ k)
    have himg2 : (e2.symm : MvPolynomial Unit (MvPolynomial n B) →ₐ[B] MvPolynomial (Unit ⊕ n) B)
        (MvPolynomial.C p * MvPolynomial.X () - 1) =
        MvPolynomial.rename Sum.inr p * MvPolynomial.X (Sum.inl ()) - 1 := by
      show e2.symm (MvPolynomial.C p * MvPolynomial.X () - 1) = _
      rw [map_sub, map_mul, map_one, he2symC, he2symX]
    rw [himg1, himg2]
  exact ⟨e0.trans (e1.trans (Ideal.quotientEquivAlg _ _ e2.symm hEq))⟩

open scoped Classical in
/-- **`flat_algEquiv`を「イデアルが別の形で与えられている」場合へ橋渡し**
——呼び出し側(`ExtLimit.lean`)では局所化の底のイデアルが
`Ideal.map (MvPolynomial.map φ) I₀`という形で与えられており、
インスタンス(`Algebra`・`IsLocalization.Away`・`IsScalarTower`)も
その形で付いている。`rw`はゴールしか書き換えないためインスタンス側と
食い違ってしまうので、イデアルを変数`Iq`として受け取り、
`Iq = Ideal.span (Set.range q₁)`という等式を仮定して`subst`する形に
しておく——こうすると呼び出し側は等式を1本渡すだけで済む。

★**sorry 無し**。標準3公理のみ。 -/
theorem localization_away_quotient_mvPolynomial_flat_algEquiv_of_eq (B' : Type) [CommRing B']
    (ι : Type) [Fintype ι] (κ₀ : Type) [Fintype κ₀]
    (Iq : Ideal (MvPolynomial ι B')) (q₁ : κ₀ → MvPolynomial ι B')
    (hIq : Iq = Ideal.span (Set.range q₁)) (p : MvPolynomial ι B')
    (S : Type) [CommRing S] [Algebra (MvPolynomial ι B' ⧸ Iq) S]
    [IsLocalization.Away (Ideal.Quotient.mk Iq p) S]
    [Algebra B' S] [IsScalarTower B' (MvPolynomial ι B' ⧸ Iq) S] :
    Nonempty (S ≃ₐ[B'] MvPolynomial (Unit ⊕ ι) B' ⧸ Ideal.span (Set.range
        (Sum.elim (fun k => MvPolynomial.rename Sum.inr (q₁ k))
          (fun _ : Unit => MvPolynomial.rename Sum.inr p * MvPolynomial.X (Sum.inl ()) - 1)))) := by
  subst hIq
  exact localization_away_quotient_mvPolynomial_flat_algEquiv B' ι κ₀ q₁ p S

open scoped Classical in
/-- **平坦化した商の中では、局所化の分母`p`の像が単元である**——
`localization_away_quotient_mvPolynomial_flat_algEquiv`が作る関係式の族に
`rename Sum.inr p * X (Sum.inl ()) - 1`が入っているので、`X (Sum.inl ())`の
クラスがそのまま逆元になる。

`t`の構成では「`D(f)`側の平坦化表示`F`が、実は`A_f[1/g]`(=`Γ(X.left,
D(f*g))`)の上の代数でもある」ことを言うのに使う——`F`の底環`A_f⊗R'`の中で
`g`に対応する元の像が単元だと分かれば、局所化の普遍性で`A_f[1/g]⊗R'`から
`F`への代数構造が誘導され、`D(f)`側と`D(g)`側が**共通の底**の上に乗る
(`exists_mvPolynomial_quotient_specIso_descend`は単一の`A`を要求するので
これが要る)。

★**sorry 無し**。標準3公理のみ。 -/
theorem isUnit_rename_of_flat_relation (B : Type) [CommRing B] (n : Type) [Fintype n]
    (κ₀ : Type) [Fintype κ₀] (q₀ : κ₀ → MvPolynomial n B) (p : MvPolynomial n B) :
    IsUnit (Ideal.Quotient.mk (Ideal.span (Set.range
        (Sum.elim (fun k => MvPolynomial.rename Sum.inr (q₀ k))
          (fun _ : Unit => MvPolynomial.rename Sum.inr p * MvPolynomial.X (Sum.inl ()) - 1))))
      (MvPolynomial.rename Sum.inr p)) := by
  set J := Ideal.span (Set.range
      (Sum.elim (fun k => MvPolynomial.rename Sum.inr (q₀ k))
        (fun _ : Unit => MvPolynomial.rename Sum.inr p * MvPolynomial.X (Sum.inl ()) - 1))) with hJ
  refine IsUnit.of_mul_eq_one (Ideal.Quotient.mk J (MvPolynomial.X (Sum.inl ()))) ?_
  rw [← map_mul, ← sub_eq_zero, ← map_one (Ideal.Quotient.mk J), ← map_sub,
    Ideal.Quotient.eq_zero_iff_mem, hJ]
  apply Ideal.subset_span
  exact ⟨Sum.inr (), rfl⟩

open scoped TensorProduct in
/-- **`exists_mvPolynomial_quotient_ringHom_descend2`の「底環が動く版」**——
出発側の底が`A`、行き先側の底が`A'`で、その間に`ℚ`-代数写像`φ : A →ₐ[ℚ] A'`が
あるときの、環準同型の降下。

`Lemma 4.1`の`GlueData`を「`V (i,j) := U_i ⊓ U_j`の片」で組む設計
(`corrhyp-goal.md`の`2026-09-05夜、続き19`)では、`f i j : V (i,j) ⟶ U i`は
「`Γ(X.left,U_i)`を底とする表示から`Γ(X.left,U_i ⊓ U_j)`を底とする表示への
**写像**」であり、底環が`φ`(制限写像)に沿って動く。既存の`descend2`は
単一の底を仮定しているが、**関係式を先に`φ`で押し出しておけば無料で
帰着できる**——それがこの補題である。

鍵の恒等式は
`(map (id A') (val R)) ∘ (map φ (id R)) = map φ (val R)`
(`Algebra.TensorProduct.map_comp`+`AlgHom.id_comp`/`comp_id`)であり、
`inclusion`版も同型に従う。

★**sorry 無し**。標準3公理のみ。 -/
theorem exists_mvPolynomial_quotient_ringHom_descend2_of_map
    (A A' : Type) [CommRing A] [Algebra ℚ A] [CommRing A'] [Algebra ℚ A'] (φ : A →ₐ[ℚ] A')
    (R R₂ : FgSubalgebra ℚ ℝ) {ι κ ι' κ' : Type} [Fintype ι] [Fintype κ] [Fintype κ']
    (q : κ → MvPolynomial ι (A ⊗[ℚ] R.1)) (q₂ : κ' → MvPolynomial ι' (A' ⊗[ℚ] R₂.1))
    (ψ : ι → MvPolynomial ι' (A' ⊗[ℚ] ℝ))
    (hψ : ∀ k, MvPolynomial.aeval ψ
        (MvPolynomial.map (Algebra.TensorProduct.map φ (Subalgebra.val R.1)).toRingHom (q k)) ∈
      Ideal.span (Set.range (fun k' => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A') (Subalgebra.val R₂.1)).toRingHom (q₂ k')))) :
    ∃ (R' : FgSubalgebra ℚ ℝ) (hR : R ≤ R') (hR₂ : R₂ ≤ R') (ev : ι → MvPolynomial ι' (A' ⊗[ℚ] R'.1)),
      (∀ i, MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A') (Subalgebra.val R'.1)).toRingHom (ev i) = ψ i) ∧
      (∀ k, MvPolynomial.aeval ev (MvPolynomial.map
          (Algebra.TensorProduct.map φ (Subalgebra.inclusion hR)).toRingHom (q k)) ∈
        Ideal.span (Set.range (fun k' => MvPolynomial.map
          (Algebra.TensorProduct.map (AlgHom.id ℚ A') (Subalgebra.inclusion hR₂)).toRingHom (q₂ k')))) := by
  set qt : κ → MvPolynomial ι (A' ⊗[ℚ] R.1) := fun k =>
    MvPolynomial.map (Algebra.TensorProduct.map φ (AlgHom.id ℚ R.1)).toRingHom (q k) with hqt
  have hcompval : ∀ (S : FgSubalgebra ℚ ℝ) (hRS : R ≤ S) (x : MvPolynomial ι (A ⊗[ℚ] R.1)),
      MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A') (Subalgebra.inclusion hRS)).toRingHom
        (MvPolynomial.map (Algebra.TensorProduct.map φ (AlgHom.id ℚ R.1)).toRingHom x)
      = MvPolynomial.map (Algebra.TensorProduct.map φ (Subalgebra.inclusion hRS)).toRingHom x := by
    intro S hRS x
    rw [MvPolynomial.map_map]
    congr 2
    have h : (Algebra.TensorProduct.map (AlgHom.id ℚ A') (Subalgebra.inclusion hRS)).comp
        (Algebra.TensorProduct.map φ (AlgHom.id ℚ R.1))
        = Algebra.TensorProduct.map φ (Subalgebra.inclusion hRS) := by
      rw [← Algebra.TensorProduct.map_comp, AlgHom.id_comp, AlgHom.comp_id]
    exact congrArg AlgHom.toRingHom h
  have hcompval' : ∀ x : MvPolynomial ι (A ⊗[ℚ] R.1),
      MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A') (Subalgebra.val R.1)).toRingHom
        (MvPolynomial.map (Algebra.TensorProduct.map φ (AlgHom.id ℚ R.1)).toRingHom x)
      = MvPolynomial.map (Algebra.TensorProduct.map φ (Subalgebra.val R.1)).toRingHom x := by
    intro x
    rw [MvPolynomial.map_map]
    congr 2
    have h : (Algebra.TensorProduct.map (AlgHom.id ℚ A') (Subalgebra.val R.1)).comp
        (Algebra.TensorProduct.map φ (AlgHom.id ℚ R.1))
        = Algebra.TensorProduct.map φ (Subalgebra.val R.1) := by
      rw [← Algebra.TensorProduct.map_comp, AlgHom.id_comp, AlgHom.comp_id]
    exact congrArg AlgHom.toRingHom h
  have hψ' : ∀ k, MvPolynomial.aeval ψ
      (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id ℚ A') (Subalgebra.val R.1)).toRingHom (qt k)) ∈
      Ideal.span (Set.range (fun k' => MvPolynomial.map
        (Algebra.TensorProduct.map (AlgHom.id ℚ A') (Subalgebra.val R₂.1)).toRingHom (q₂ k'))) := by
    intro k
    rw [hqt, hcompval']
    exact hψ k
  obtain ⟨R', hR, hR₂, ev, hev, hmem⟩ :=
    exists_mvPolynomial_quotient_ringHom_descend2 A' R R₂ qt q₂ ψ hψ'
  refine ⟨R', hR, hR₂, ev, hev, ?_⟩
  intro k
  have hk := hmem k
  rw [hqt] at hk
  rw [hcompval R' hR (q k)] at hk
  exact hk

/-- **制限写像が代数構造と可換なら、関係式は「移した先」でも消える**——
`Lemma 4.1`の新設計で`f i j`を`exists_mvPolynomial_quotient_ringHom_descend2_of_map`
により降ろすとき、その仮説`hψ`(関係式の像がイデアルに入ること)を検証する
**計算の核**である。`CorrHyp`に依存しない純粋な多項式環の事実として書いた。

読み方: `𝔸 = A_U⊗ℝ`・`𝔹 = A_V⊗ℝ`・`S = Γ(C,piece(U))`・`T = Γ(C,piece(V))`、
`φ'`が底の写像、`algU`/`algV`が代数構造、`restr`が制限、`valU`/`valV`が
`Presentation`の生成元、`ψ`が「`valU`の像を`𝔹`係数の多項式へ持ち上げたもの」。
`hcomm`(代数構造の可換性、`ExtLimit.lean`の`pieceAlgebraMap_naturality`)と
`hψval`(`ψ`の作り方)から、`U`側で消える関係式`p`は`V`側でも消える。

証明はmathlibの`eval₂_comp_left`(環準同型を`eval₂`の中へ入れる)と
`eval₂_map`(係数写像を`eval₂`へ吸収する)を2回ずつ使うだけ。

★**sorry 無し**。標準3公理のみ。 -/
theorem eval₂_map_aeval_eq_zero {𝔸 𝔹 S T : Type} [CommRing 𝔸] [CommRing 𝔹] [CommRing S] [CommRing T]
    (φ' : 𝔸 →+* 𝔹) (algU : 𝔸 →+* S) (algV : 𝔹 →+* T) (restr : S →+* T)
    (hcomm : restr.comp algU = algV.comp φ')
    {ι ι' : Type} (valU : ι → S) (valV : ι' → T) (ψ : ι → MvPolynomial ι' 𝔹)
    (hψval : ∀ i, MvPolynomial.eval₂ algV valV (ψ i) = restr (valU i))
    (p : MvPolynomial ι 𝔸) (hp : MvPolynomial.eval₂ algU valU p = 0) :
    MvPolynomial.eval₂ algV valV
      (MvPolynomial.eval₂ (MvPolynomial.C : 𝔹 →+* MvPolynomial ι' 𝔹) ψ (MvPolynomial.map φ' p)) = 0 := by
  show (MvPolynomial.eval₂Hom algV valV)
    (MvPolynomial.eval₂ (MvPolynomial.C : 𝔹 →+* MvPolynomial ι' 𝔹) ψ (MvPolynomial.map φ' p)) = 0
  rw [MvPolynomial.eval₂_comp_left (MvPolynomial.eval₂Hom algV valV)]
  have hC : (MvPolynomial.eval₂Hom algV valV).comp (MvPolynomial.C : 𝔹 →+* MvPolynomial ι' 𝔹) = algV := by
    ext b
    simp
  have hfun : ((MvPolynomial.eval₂Hom algV valV : MvPolynomial ι' 𝔹 →+* T) ∘ ψ)
      = (⇑restr ∘ valU) := funext hψval
  rw [hC, hfun, MvPolynomial.eval₂_map, ← hcomm,
    ← MvPolynomial.eval₂_comp_left restr algU valU p, hp, map_zero]

/-- **`= 0` から `∈ J` へ戻す**——表示(presentation)があれば
`eval₂`の核はちょうど関係式のイデアルである。

`e : MvPolynomial ι' 𝔹 ⧸ J ≃+* T` が「`T` は生成元 `X i` と係数 `C b` で
`J` を法として表示される」という情報を持っているとき、`eval₂ algV valV` は
`e ∘ mk` と**環準同型として一致する**(`ringHom_ext`で生成元と定数だけ見れば
よい)。よって `eval₂ algV valV r = 0` は `mk r = 0`、すなわち `r ∈ J`。

`eval₂_map_aeval_eq_zero`の相棒。あちらが「消える」ことを言い、
こちらが「消えるならイデアルに入る」ことを言う。

★**sorry 無し**。標準3公理のみ。 -/
theorem mem_ideal_of_eval₂_eq_zero {𝔹 T : Type} [CommRing 𝔹] [CommRing T] {ι' : Type}
    (J : Ideal (MvPolynomial ι' 𝔹)) (e : (MvPolynomial ι' 𝔹 ⧸ J) ≃+* T)
    (algV : 𝔹 →+* T) (valV : ι' → T)
    (halg : ∀ b, algV b = e (Ideal.Quotient.mk J (MvPolynomial.C b)))
    (hval : ∀ i, valV i = e (Ideal.Quotient.mk J (MvPolynomial.X i)))
    (r : MvPolynomial ι' 𝔹) (hr : MvPolynomial.eval₂ algV valV r = 0) : r ∈ J := by
  have hEq : (MvPolynomial.eval₂Hom algV valV)
      = ((e : (MvPolynomial ι' 𝔹 ⧸ J) →+* T).comp (Ideal.Quotient.mk J)) := by
    apply MvPolynomial.ringHom_ext
    · intro b; simpa using halg b
    · intro i; simpa using hval i
  have h0 : e (Ideal.Quotient.mk J r) = 0 := by
    have hh : (MvPolynomial.eval₂Hom algV valV) r
        = ((e : (MvPolynomial ι' 𝔹 ⧸ J) →+* T).comp (Ideal.Quotient.mk J)) r :=
      congrArg (fun f : MvPolynomial ι' 𝔹 →+* T => f r) hEq
    simpa [hr] using hh.symm
  rw [← Ideal.Quotient.eq_zero_iff_mem]
  exact e.map_eq_zero_iff.mp h0

/-- **`hψ`そのもの**——2本を合成した、`descend2_of_map`にそのまま食わせる形。

`U`側で消える関係式 `p` を、係数を `φ'` で `V` 側へ移し生成元に `ψ` を代入
すると、`V`側の表示のイデアル `J` に入る。`Lemma 4.1`の`f i j`を降ろすとき、
`hcomm`には`ExtLimit.lean`の`pieceAlgebraMap_naturality`が、`J`には
`exists_descendPieceR_flat_mvPolynomial_baseChange`が与える平坦表示が入る。

★**sorry 無し**。標準3公理のみ。 -/
theorem aeval_map_mem_ideal_of_relation
    {𝔸 𝔹 S T : Type} [CommRing 𝔸] [CommRing 𝔹] [CommRing S] [CommRing T]
    {ι ι' : Type}
    (φ' : 𝔸 →+* 𝔹) (algU : 𝔸 →+* S) (algV : 𝔹 →+* T) (restr : S →+* T)
    (hcomm : restr.comp algU = algV.comp φ')
    (valU : ι → S) (valV : ι' → T) (ψ : ι → MvPolynomial ι' 𝔹)
    (hψval : ∀ i, MvPolynomial.eval₂ algV valV (ψ i) = restr (valU i))
    (J : Ideal (MvPolynomial ι' 𝔹)) (e : (MvPolynomial ι' 𝔹 ⧸ J) ≃+* T)
    (halg : ∀ b, algV b = e (Ideal.Quotient.mk J (MvPolynomial.C b)))
    (hval : ∀ i, valV i = e (Ideal.Quotient.mk J (MvPolynomial.X i)))
    (p : MvPolynomial ι 𝔸) (hp : MvPolynomial.eval₂ algU valU p = 0) :
    MvPolynomial.aeval ψ (MvPolynomial.map φ' p) ∈ J := by
  rw [MvPolynomial.aeval_def, MvPolynomial.algebraMap_eq]
  exact mem_ideal_of_eval₂_eq_zero J e algV valV halg hval _
    (eval₂_map_aeval_eq_zero φ' algU algV restr hcomm valU valV ψ hψval p hp)

/-- **`Algebra.Generators` があれば `ψ` は自分で作れる**——
`Generators` は生成元 `val` だけでなく**切断 `σ`** (`aeval val (σ s) = s`)
を持つので、`ψ i := PV.σ (restr (PU.val i))` と置けば
`eval₂_map_aeval_eq_zero` の `hψval` は `Generators.aeval_val_σ` そのもの。
つまり**`hcomm`(代数構造が制限と可換)さえあれば `hψ` は自動**である。

`Lemma 4.1` の `f i j` を降ろすとき、`hcomm` には `ExtLimit.lean` の
`pieceAlgebraMap_naturality` がそのまま入る。

★**sorry 無し**。標準3公理のみ。 -/
theorem aeval_map_relation_mem_ker
    {𝔸 𝔹 S T : Type} [CommRing 𝔸] [CommRing 𝔹] [CommRing S] [CommRing T]
    [Algebra 𝔸 S] [Algebra 𝔹 T]
    {ι ι' : Type} (PU : Algebra.Generators 𝔸 S ι) (PV : Algebra.Generators 𝔹 T ι')
    (φ' : 𝔸 →+* 𝔹) (restr : S →+* T)
    (hcomm : restr.comp (algebraMap 𝔸 S) = (algebraMap 𝔹 T).comp φ')
    (p : MvPolynomial ι 𝔸) (hp : p ∈ PU.ker) :
    MvPolynomial.aeval (fun i => PV.σ (restr (PU.val i))) (MvPolynomial.map φ' p) ∈ PV.ker := by
  rw [Algebra.Generators.ker_eq_ker_aeval_val, RingHom.mem_ker]
  simp only [MvPolynomial.aeval_def, MvPolynomial.algebraMap_eq]
  have hp0 : MvPolynomial.eval₂ (algebraMap 𝔸 S) PU.val p = 0 := by
    have h := hp
    rw [Algebra.Generators.ker_eq_ker_aeval_val, RingHom.mem_ker] at h
    simpa [MvPolynomial.aeval_def] using h
  have hψval : ∀ i, MvPolynomial.eval₂ (algebraMap 𝔹 T) PV.val (PV.σ (restr (PU.val i)))
      = restr (PU.val i) := by
    intro i
    have h := PV.aeval_val_σ (restr (PU.val i))
    simpa [MvPolynomial.aeval_def] using h
  exact eval₂_map_aeval_eq_zero φ' (algebraMap 𝔸 S) (algebraMap 𝔹 T) restr hcomm
    PU.val PV.val (fun i => PV.σ (restr (PU.val i))) hψval p hp0

/-- **`descend2_of_map` の `hψ` の最終形**——`Presentation` を2つ与えれば、
「`U` 側の関係式 `k` を `V` 側へ移したものは `V` 側の関係式が張るイデアルに
入る」が `hcomm` だけから出る。`Presentation` の構造フィールド
`span_range_relation_eq_ker` で `ker` と `span (range relation)` を
往復するだけで `aeval_map_relation_mem_ker` に帰着する。

★**sorry 無し**。標準3公理のみ。 -/
theorem aeval_map_relation_mem_span
    {𝔸 𝔹 S T : Type} [CommRing 𝔸] [CommRing 𝔹] [CommRing S] [CommRing T]
    [Algebra 𝔸 S] [Algebra 𝔹 T]
    {ι κ ι' κ' : Type} (PU : Algebra.Presentation 𝔸 S ι κ) (PV : Algebra.Presentation 𝔹 T ι' κ')
    (φ' : 𝔸 →+* 𝔹) (restr : S →+* T)
    (hcomm : restr.comp (algebraMap 𝔸 S) = (algebraMap 𝔹 T).comp φ') (k : κ) :
    MvPolynomial.aeval (fun i => PV.σ (restr (PU.val i)))
        (MvPolynomial.map φ' (PU.relation k))
      ∈ Ideal.span (Set.range PV.relation) := by
  rw [PV.span_range_relation_eq_ker]
  refine aeval_map_relation_mem_ker PU.toGenerators PV.toGenerators φ' restr hcomm _ ?_
  rw [← PU.span_range_relation_eq_ker]
  exact Ideal.subset_span ⟨k, rfl⟩

open scoped TensorProduct in
/-- **`f ⊗ g` を「係数だけ動かす」と「底だけ動かす」に分ける**——
`Lemma 4.1`で`R`レベルの関係式`q₀`(`A ⊗ R.1`係数)を`A' ⊗ ℝ`へ移すとき、
`φ ⊗ (val R)` を `(φ ⊗ id) ∘ (id ⊗ val R)` と分解すると、後半が
`pieceAlgebra_relation_descend_q₀_map`(降ろした関係式の定義方程式)、
前半が`aeval_map_relation_mem_span`の`φ'`にそのまま対応する。

★**sorry 無し**。標準3公理のみ。 -/
theorem tensorProduct_map_split {R A A' B B' : Type} [CommRing R] [CommRing A] [CommRing A']
    [CommRing B] [CommRing B'] [Algebra R A] [Algebra R A'] [Algebra R B] [Algebra R B']
    (f : A →ₐ[R] A') (g : B →ₐ[R] B') :
    Algebra.TensorProduct.map f g
      = (Algebra.TensorProduct.map f (AlgHom.id R B')).comp
          (Algebra.TensorProduct.map (AlgHom.id R A) g) := by
  ext <;> simp

set_option maxHeartbeats 2000000 in
open scoped TensorProduct in
/-- `tensorProduct_map_split`の`MvPolynomial`版。**使う側でこの形が要る**
ので、`CorrHyp`の巨大な型を持ち込まない小さい文脈でここに置いておく
(`Algebra.TensorProduct`のインスタンス解決が重く、ここで56秒かかるが、
使う側では`exact`一発・3秒で済む)。

★**sorry 無し**。標準3公理のみ。 -/
theorem mvPolynomial_map_tensorMap_split {R A A' B B' : Type} [CommRing R] [CommRing A]
    [CommRing A'] [CommRing B] [CommRing B'] [Algebra R A] [Algebra R A'] [Algebra R B]
    [Algebra R B'] (f : A →ₐ[R] A') (g : B →ₐ[R] B') {ι : Type} (x : MvPolynomial ι (A ⊗[R] B)) :
    MvPolynomial.map (Algebra.TensorProduct.map f (AlgHom.id R B')).toRingHom
      (MvPolynomial.map (Algebra.TensorProduct.map (AlgHom.id R A) g).toRingHom x)
    = MvPolynomial.map (Algebra.TensorProduct.map f g).toRingHom x :=
  Eq.trans (MvPolynomial.map_map _ _ x)
    (congrArg (fun t : (A ⊗[R] B) →ₐ[R] (A' ⊗[R] B') => MvPolynomial.map t.toRingHom x)
      (tensorProduct_map_split f g).symm)

/-- **商環の間の環準同型を、生成元の行き先`ev`と係数の写像`ρ`から作る**——
`Lemma 4.1`の`f i j`は「`R'`レベルの環準同型の`Spec`」なので、
`descend2_of_map`が返す生データ(`ev`と関係式の所属)から実際の`RingHom`を
組み立てるのがこれ。`Ideal.Quotient.lift`のwell-definednessは
`Ideal.span_le`＋`Ideal.comap`で`hq`から出る。

既存の`exists_mvPolynomial_quotient_ringEquiv_of_data`は**同型**を作る
ものだったが、新設計(`V (i,j) := U_i ⊓ U_j`の片)では`f i j`は**写像**で
よく、係数環も`T`から`T'`へ動く。その分だけ軽い。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def mvPolynomial_quotient_ringHom_of_data {T T' : Type} [CommRing T] [CommRing T']
    {ι κ ι' κ' : Type} (ρ : T →+* T')
    (q : κ → MvPolynomial ι T) (q₂ : κ' → MvPolynomial ι' T') (ev : ι → MvPolynomial ι' T')
    (hq : ∀ k, MvPolynomial.aeval ev (MvPolynomial.map ρ (q k)) ∈ Ideal.span (Set.range q₂)) :
    (MvPolynomial ι T ⧸ Ideal.span (Set.range q)) →+*
      (MvPolynomial ι' T' ⧸ Ideal.span (Set.range q₂)) :=
  Ideal.Quotient.lift _
    ((Ideal.Quotient.mk (Ideal.span (Set.range q₂))).comp
      ((MvPolynomial.aeval ev).toRingHom.comp (MvPolynomial.map ρ)))
    (fun a ha => by
      rw [RingHom.comp_apply, Ideal.Quotient.eq_zero_iff_mem]
      have hle : Ideal.span (Set.range q) ≤ Ideal.comap
          ((MvPolynomial.aeval ev).toRingHom.comp (MvPolynomial.map ρ))
          (Ideal.span (Set.range q₂)) :=
        Ideal.span_le.mpr (fun x hx => by obtain ⟨k, rfl⟩ := hx; exact hq k)
      exact hle ha)

/-- `mvPolynomial_quotient_ringHom_of_data`の計算則(`rfl`)。 -/
theorem mvPolynomial_quotient_ringHom_of_data_mk {T T' : Type} [CommRing T] [CommRing T']
    {ι κ ι' κ' : Type} (ρ : T →+* T')
    (q : κ → MvPolynomial ι T) (q₂ : κ' → MvPolynomial ι' T') (ev : ι → MvPolynomial ι' T')
    (hq : ∀ k, MvPolynomial.aeval ev (MvPolynomial.map ρ (q k)) ∈ Ideal.span (Set.range q₂))
    (p : MvPolynomial ι T) :
    mvPolynomial_quotient_ringHom_of_data ρ q q₂ ev hq
        (Ideal.Quotient.mk (Ideal.span (Set.range q)) p)
      = Ideal.Quotient.mk (Ideal.span (Set.range q₂))
        (MvPolynomial.aeval ev (MvPolynomial.map ρ p)) := rfl

end ABC3.Found.CorrHyp
