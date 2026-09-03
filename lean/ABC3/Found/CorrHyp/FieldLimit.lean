/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.RingTheory.Adjoin.FG
import Mathlib.Algebra.Category.Ring.Basic
import Mathlib.CategoryTheory.Category.Preorder
import Mathlib.AlgebraicGeometry.Scheme
import Mathlib.AlgebraicGeometry.GammaSpecAdjunction
import Mathlib.AlgebraicGeometry.AffineTransitionLimit

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

/- ★★次の一手(未着手): `Lemma 4.1` 本体へ——`X_K`・`Z_K`・その間の
correspondence を「`toSchemeDiagram` 上の適当な `i`(有限生成部分環 `R`)への
下降」として表すには、`Scheme.exists_hom_comp_eq_comp_of_locallyOfFiniteType`
等の側条件(`∀ {i j} (f : i ⟶ j), IsAffineHom (D.map f)`——`Spec S → Spec R`
は常にアフィン射なので無条件で成り立つはず、`[∀ i, CompactSpace/QuasiSeparatedSpace
(D.obj i)]`——`Spec R` は常にコンパクト・準分離、`LocallyOfFinitePresentation f`
——`X_K → Spec K` 側に要る仮定)を揃え、`HyperbolicCurveData` の
`IsGenericallyScheme`/`ModuliStack` 関連フィールドと接続する。 -/

end ABC3.Found.CorrHyp
