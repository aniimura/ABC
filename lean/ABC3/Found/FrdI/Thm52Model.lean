/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm52BiratSec
import ABC3.Found.FrdI.Cor411Otimes
import ABC3.Found.FrdI.Prop44Phi
import ABC3.Found.FrdI.Thm52

/-!
# [FrdI] Theorem 5.2, (iv) の第 4 段 —— 関手 `𝒞̃ ⥤ 𝒞^model`

原文 (FrdI p.102):
> verify that these assignments determine a functor

★★本ファイルで **`Theorem 5.2, (iv)` の関手そのもの**を作る。

* 対象 `X ↦ (Base X, cls X)` —— `cls` は `F-𝒫-path` から決まる(`Thm52Path.lean`)。
* 射 `f ↦ (deg f, Base f, Div f, u_f)` —— `u_f` が (iv) の**核心**。

★★`u_f` の作り方(原文どおり):
1. `w_f := π_X⁻¹ ≫ f^birat ≫ π_Y` を取る(`Thm52BiratSec.lean`)。
2. `Remark 2.7.2`(等長版)で `w_f = F_n ≫ β_f ≫ a_f` と一意分解する。
3. `u_f := B(Base π_X)(κ(β_f))`。

★★条件式が合うのは `pathW_divGp`(`w_f` の零因子)がちょうど
`deg(f)·cls X + Div(f) − Φ(Base f)(cls Y)` になるからであり、
★★合成則が合うのは `rem272Beta_comp`(中央項の合成則)によるものである。

★**有理関数の単系 `B`** は原文が (iv) の仮定として置くもの
(「`B` は `𝒞` に伴う有理関数の単系」)なので、Lean でも
**interface(`RatFnData`)として仮定に置く**。必要なのは
`B ≅ 𝒪^×(−^birat)`(`kappa`)、`Div_B ≅ Div^gp`(`kappa_divB`)、
`𝒟` の射に沿った自然性(`kappa_pull`)の 3 点だけである。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-! ## ★1. `path` が定める底の同型 -/

/-- ★★`F-𝒫-path` が定める **底の同型** `Base X ≅ Base ref`。

★`biratBase (piHom …)` は `biratDown` を経由した型を持ち、`𝒞` 側の型と
**定義上は等しいが構文上は違う**ため、`rw`/`simp` が通らない。
★そこで **`𝒞` 側の型で書いた別名**を置く。 -/
noncomputable def FPPath.piBase {S : BaseSection P} {X : C} (p : FPPath S X) :
    (P.toElem.obj X).base ⟶ (P.toElem.obj p.ref).base :=
  @inv _ _ _ _ (P.Base p.toObj) p.toObj_preStep.2 ≫ P.Base p.toRef

instance FPPath.isIso_piBase {S : BaseSection P} {X : C} (p : FPPath S X) :
    IsIso p.piBase := by
  haveI h1 : IsIso (P.Base p.toObj) := p.toObj_preStep.2
  haveI h2 : IsIso (P.Base p.toRef) := p.toRef_preStep.2
  show IsIso (@inv _ _ _ _ (P.Base p.toObj) p.toObj_preStep.2 ≫ P.Base p.toRef)
  infer_instance

theorem FPPath.biratBase_piHom (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {S : BaseSection P} {X : C} (p : FPPath S X) :
    biratBase (p.piHom G hiso) = p.piBase := p.biratBase_biratIso G hiso

/-- ★★`w_f` の底を `𝒟` の中で書いたもの。 -/
noncomputable def pathTheta {S : BaseSection P} {X Y : PathCat S} (f : X ⟶ Y) :
    (P.toElem.obj X.path.ref).base ⟶ (P.toElem.obj Y.path.ref).base :=
  inv X.path.piBase ≫ P.Base f ≫ Y.path.piBase

/-- ★★★**`u` の合成則の骨格** —— `π_X ≫ θ_f = Base f ≫ π_Y`。 -/
theorem pathTheta_spec {S : BaseSection P} {X Y : PathCat S} (f : X ⟶ Y) :
    X.path.piBase ≫ pathTheta f = P.Base f ≫ Y.path.piBase := by
  rw [pathTheta, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp]

/-- ★★`w_f` の底は `pathTheta`。 -/
theorem biratBase_pathW (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {S : BaseSection P} {X Y : PathCat S} (f : X ⟶ Y) :
    biratBase (pathW G hiso f) = pathTheta f := by
  have hX : biratBase (X.path.piHom G hiso) = X.path.piBase :=
    FPPath.biratBase_piHom G hiso X.path
  have hY : biratBase (Y.path.piHom G hiso) = Y.path.piBase :=
    FPPath.biratBase_piHom G hiso Y.path
  have h : X.path.piBase ≫ biratBase (pathW G hiso f) = P.Base f ≫ Y.path.piBase := by
    refine Eq.trans ?_ (Eq.trans (pathW_base G hiso f) ?_)
    · exact congrArg (fun m : (P.toElem.obj X.obj).base ⟶ (P.toElem.obj X.path.ref).base =>
        m ≫ biratBase (pathW G hiso f)) hX.symm
    · exact congrArg (fun m : (P.toElem.obj Y.obj).base ⟶ (P.toElem.obj Y.path.ref).base =>
        P.Base f ≫ m) hY
  rw [pathTheta, ← h]
  exact (IsIso.inv_hom_id_assoc X.path.piBase _).symm

/-- ★`pathW_divGp` を `𝒞` 側の型で書き直したもの。 -/
theorem pathW_divGp' (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {S : BaseSection P} {X Y : PathCat S} (f : X ⟶ Y) :
    Φ.gpMapOn X.path.piBase (biratDivGp (pathW G hiso f))
      = ((P.degFr f : ℕ+) : ℕ) • X.path.cls + toGp _ (P.Div f)
        - Φ.gpMapOn (P.Base f) Y.path.cls := by
  refine Eq.trans ?_ (pathW_divGp G hiso f)
  exact congrArg (fun m : (P.toElem.obj X.obj).base ⟶ (P.toElem.obj X.path.ref).base =>
    Φ.gpMapOn m (biratDivGp (pathW G hiso f))) (FPPath.biratBase_piHom G hiso X.path).symm

/-! ## ★2. 有理関数の単系 `B` の interface

原文 (FrdI p.101):
> (ii) The category C is a Frobenioid [with respect to the functor C → FΦ] of isotropic and

★原文は (iv) の仮定として「`B` は `𝒞` に伴う有理関数の単系」と置き、
(ii) で「`𝒪^×(−)` と `B` の間に自然同型があり、`Div_B` と両立する」と言う。
★**そのまま構造として書く**のが忠実である。 -/

/-- ★★**有理関数の単系 `B` の interface**。

* `kappa` —— `𝒪^×(A^birat) ≅ B(Base A)`(原文 (ii) の「自然同型」)
* `kappa_divB` —— `Div_B` と `Div^gp` の両立(原文 (ii) の「compatible」)
* `kappa_pull` —— `𝒟` の射(pull-back 射の底)に沿った自然性 -/
structure RatFnData (P : PreFrobenioid C Φ) (G : Frobenioid P) where
  /-- 有理関数の単系 `B` -/
  bmon : MonoidOn.{v, u, w} D
  /-- `Div_B : B → Φ^gp` -/
  divB : ∀ A : D, bmon.val A →+ Gp (Φ.val A)
  /-- `Div_B` の自然性 -/
  divB_nat : ∀ {A B : D} (f : A ⟶ B) (x : bmon.val B),
    divB A (bmon.map f x) = Φ.gpMapOn f (divB B x)
  /-- `𝒪^×(A^birat) ≅ B(Base A)` -/
  kappa : ∀ A : C, Additive ↥(OTimes (biratPre P G) ((toBiratCat P G).obj A))
    ≃+ bmon.val (P.toElem.obj A).base
  /-- `Div_B ∘ κ = Div^gp` -/
  kappa_divB : ∀ (A : C) (δ : OTimes (biratPre P G) ((toBiratCat P G).obj A)),
    divB _ (kappa A (Additive.ofMul δ))
      = biratDivGp ((δ : End ((toBiratCat P G).obj A)) : _ ⟶ _)
  /-- pull-back 射に沿った自然性 -/
  kappa_pull : ∀ {A B : C} (φ : (toBiratCat P G).obj A ⟶ (toBiratCat P G).obj B)
    (θ : (P.toElem.obj A).base ⟶ (P.toElem.obj B).base)
    (δ : OTimes (biratPre P G) ((toBiratCat P G).obj B))
    (ε : OTimes (biratPre P G) ((toBiratCat P G).obj A)),
    IsPullBack (biratPre P G) φ → biratBase φ = θ →
    φ ≫ ((δ : End ((toBiratCat P G).obj B)) : _ ⟶ _)
      = ((ε : End ((toBiratCat P G).obj A)) : _ ⟶ _) ≫ φ →
    kappa A (Additive.ofMul ε) = bmon.map θ (kappa B (Additive.ofMul δ))

/-- ★interface から `ModelData` を作る(`Φ` はそのまま)。 -/
def RatFnData.model {G : Frobenioid P} (R : RatFnData P G) : ModelData.{v, u, w} D where
  phi := Φ
  bmon := R.bmon
  divB := R.divB
  divB_nat := R.divB_nat

/-- ★`𝒞^birat` の Frobenioid 核。 -/
noncomputable def biratCoreOf (G : Frobenioid P)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X) :
    FrobenioidCore (biratPre P G) := (biratFrobenioid P G hfn).core

/-- ★`𝒞^birat` はすべて等長。 -/
theorem biratMet (G : Frobenioid P) {A B : BiratCat P G} (f : A ⟶ B) :
    IsIsometric (biratPre P G) f := birat_isIsometric f

/-! ## ★3. 単元成分 `u_f` -/

variable [IsConnected D]

/-- ★★**`w_f` の 3 分解の中央項** `β_f ∈ 𝒪^×(ref_X^birat)`。 -/
noncomputable def pathBetaO (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs)
    {X Y : PathCat S} (f : X ⟶ Y) :
    OTimes (biratPre P G) ((toBiratCat P G).obj X.path.ref) :=
  rem272BetaO (biratCoreOf G hfn) (fun Z => birat_isOfIsotropicType hiso Z)
    (fun {_ _} q => biratMet G q) (biratFrobSection_isFrobeniusSection G hiso S hFs)
    X.path.ref_mem Y.path.ref_mem (pathW G hiso f)

/-- ★★★**model 射の単元成分 `u_f`**。

原文 (FrdI p.102):
> of the respective objects and morphisms of C in Cbirat] of C. Now in light of the fact
-/
noncomputable def pathU {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs)
    {X Y : PathCat S} (f : X ⟶ Y) : R.bmon.val (P.toElem.obj X.obj).base :=
  R.bmon.map X.path.piBase
    (R.kappa X.path.ref (Additive.ofMul (pathBetaO G hiso hfn S Fs hFs f)))

theorem pathU_eq {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs)
    {X Y : PathCat S} (f : X ⟶ Y) :
    pathU R hiso hfn S Fs hFs f
      = R.bmon.map X.path.piBase
          (R.kappa X.path.ref (Additive.ofMul (pathBetaO G hiso hfn S Fs hFs f))) := rfl

/-! ## ★4. `Div_B(u_f)` の計算 -/

/-- ★`𝒫^birat` の射の `Div^gp` は 0。 -/
theorem biratDivGp_of_homP (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (S : BaseSection P) {A B : BiratCat P G} {g : A ⟶ B}
    (hg : (biratBaseSection G hiso S).homP g) : biratDivGp g = 0 := by
  obtain ⟨f₀, hf, rfl⟩ := hg
  refine (biratDivGp_toBiratMap G f₀).trans ?_
  rw [show P.Div f₀ = 0 from (G.core.pullBackLB f₀ (S.isPullBack hf)).1.2]
  exact toGp_zero _

/-- ★`𝒞^birat` の Frobenius-section の `Div^gp` は 0。 -/
theorem biratDivGp_frobSection (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (S : BaseSection P) {Fs : ℕ+ →* SectionEnd S} (hFs : IsFrobeniusSection S Fs)
    (n : ℕ+) (A : (biratBaseSection G hiso S).Obj) :
    biratDivGp ((((biratFrobSection G hiso S Fs) n).app A : End A.1) : A.1 ⟶ A.1) = 0 := by
  refine (biratDivGp_toBiratMap G
    (((Fs n).app ⟨biratDown P G A.1, A.2⟩ : End _) : _ ⟶ _)).trans ?_
  rw [show P.Div (((Fs n).app ⟨biratDown P G A.1, A.2⟩ : End _) : _ ⟶ _)
    = 0 from (hFs.frobType n ⟨_, A.2⟩).1.2]
  exact toGp_zero _

/-- ★★★**3 分解の中央項の `Div^gp` は、もとの射の `Div^gp` に等しい**
—— `F` と `a` の `Div^gp` が 0 だから。 -/
theorem biratDivGp_rem272Beta (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs)
    {A B : BiratCat P G} (hA : (biratBaseSection G hiso S).objP A)
    (hB : (biratBaseSection G hiso S).objP B) (φ : A ⟶ B) :
    biratDivGp (rem272Beta (biratCoreOf G hfn) (fun Z => birat_isOfIsotropicType hiso Z)
        (fun {_ _} q => biratMet G q) (biratFrobSection_isFrobeniusSection G hiso S hFs)
        hA hB φ) = biratDivGp φ := by
  set S' := biratBaseSection G hiso S with hS'
  set Fs' := biratFrobSection G hiso S Fs with hFs'def
  have hFs'' := biratFrobSection_isFrobeniusSection G hiso S hFs
  set β := rem272Beta (biratCoreOf G hfn) (fun Z => birat_isOfIsotropicType hiso Z)
      (fun {_ _} q => biratMet G q) hFs'' hA hB φ with hβdef
  have hspec := rem272Beta_spec (biratCoreOf G hfn) (fun Z => birat_isOfIsotropicType hiso Z)
      (fun {_ _} q => biratMet G q) hFs'' hA hB φ
  set n := (biratPre P G).degFr φ with hn
  set Fn := (((Fs' n).app ⟨A, hA⟩ : End A) : A ⟶ A) with hFn
  set α := S'.lift hA hB ((biratPre P G).Base φ) with hα
  have hdα : biratDivGp α = 0 := biratDivGp_of_homP G hiso S (S'.lift_homP hA hB _)
  have hdFn : biratDivGp Fn = 0 := biratDivGp_frobSection G hiso S hFs n ⟨A, hA⟩
  have hdegα : biratDeg α = 1 :=
    ((biratCoreOf G hfn).pullBackLB α (S'.isPullBack (S'.lift_homP hA hB _))).2
  have hbFn : IsBaseIdentity (biratPre P G) Fn := hFs''.baseIdentity n ⟨A, hA⟩
  have h := congrArg biratDivGp hspec
  rw [biratDivGp_comp', biratDivGp_comp', hdα, hdFn, hdegα, map_zero, zero_add,
    PNat.one_coe, one_smul, smul_zero, add_zero,
    gpMap_biratBase_of_baseIdentity hbFn] at h
  exact h.symm

/-- ★★★★**`Div_B(u_f)` の計算** —— これが model Frobenioid の条件式そのもの。 -/
theorem pathU_divB {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs)
    {X Y : PathCat S} (f : X ⟶ Y) :
    R.divB _ (pathU R hiso hfn S Fs hFs f)
      = ((P.degFr f : ℕ+) : ℕ) • X.path.cls + toGp _ (P.Div f)
        - Φ.gpMapOn (P.Base f) Y.path.cls := by
  have hinner : R.divB _ (R.kappa X.path.ref
        (Additive.ofMul (pathBetaO G hiso hfn S Fs hFs f)))
      = biratDivGp (pathW G hiso f) :=
    (R.kappa_divB X.path.ref (pathBetaO G hiso hfn S Fs hFs f)).trans
      (biratDivGp_rem272Beta G hiso hfn S Fs hFs X.path.ref_mem Y.path.ref_mem
        (pathW G hiso f))
  refine ((R.divB_nat X.path.piBase
      (R.kappa X.path.ref (Additive.ofMul (pathBetaO G hiso hfn S Fs hFs f)))).trans
    (congrArg (⇑(Φ.gpMapOn X.path.piBase)) hinner)).trans (pathW_divGp' G hiso f)

/-! ## ★5. `u` の合成則 -/

/-- ★`w_f` の 3 分解の `P`-distinguished 成分 `a_f`。 -/
noncomputable def pathAf (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (S : BaseSection P) {X Y : PathCat S} (f : X ⟶ Y) :
    (toBiratCat P G).obj X.path.ref ⟶ (toBiratCat P G).obj Y.path.ref :=
  (biratBaseSection G hiso S).lift X.path.ref_mem Y.path.ref_mem
    ((biratPre P G).Base (pathW G hiso f))

theorem pathAf_isPullBack (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (S : BaseSection P) {X Y : PathCat S} (f : X ⟶ Y) :
    IsPullBack (biratPre P G) (pathAf G hiso S f) :=
  (biratBaseSection G hiso S).isPullBack
    ((biratBaseSection G hiso S).lift_homP X.path.ref_mem Y.path.ref_mem _)

theorem pathAf_base (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (S : BaseSection P) {X Y : PathCat S} (f : X ⟶ Y) :
    biratBase (pathAf G hiso S f) = pathTheta f :=
  ((biratBaseSection G hiso S).lift_base X.path.ref_mem Y.path.ref_mem _).trans
    (biratBase_pathW G hiso f)

/-- ★★`γ_{f,g}` —— `β_g` を `a_f` に沿って引き戻したもの。 -/
noncomputable def pathGamma (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs)
    {X Y Z : PathCat S} (f : X ⟶ Y) (g : Y ⟶ Z) :
    OTimes (biratPre P G) ((toBiratCat P G).obj X.path.ref) :=
  ⟨otimesPull (biratPre P G) (biratCoreOf G hfn) (fun W => birat_isOfIsotropicType hiso W)
      (pathAf G hiso S f) (pathAf_isPullBack G hiso S f) (pathBetaO G hiso hfn S Fs hFs g),
    otimesPull_mem (biratPre P G) (biratCoreOf G hfn) (fun W => birat_isOfIsotropicType hiso W)
      (pathAf G hiso S f) (pathAf_isPullBack G hiso S f) (pathBetaO G hiso hfn S Fs hFs g)⟩

theorem pathGamma_spec (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs)
    {X Y Z : PathCat S} (f : X ⟶ Y) (g : Y ⟶ Z) :
    pathAf G hiso S f ≫ ((pathBetaO G hiso hfn S Fs hFs g :
        End ((toBiratCat P G).obj Y.path.ref)) : _ ⟶ _)
      = ((pathGamma G hiso hfn S Fs hFs f g :
          End ((toBiratCat P G).obj X.path.ref)) : _ ⟶ _) ≫ pathAf G hiso S f :=
  otimesPull_spec (biratPre P G) (biratCoreOf G hfn) (fun W => birat_isOfIsotropicType hiso W)
    (pathAf G hiso S f) (pathAf_isPullBack G hiso S f) (pathBetaO G hiso hfn S Fs hFs g)

/-- ★★★★**`β` の合成則**(path つきの圏で)。 -/
theorem pathBetaO_comp (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs)
    {X Y Z : PathCat S} (f : X ⟶ Y) (g : Y ⟶ Z) :
    pathBetaO G hiso hfn S Fs hFs (f ≫ g)
      = pathGamma G hiso hfn S Fs hFs f g
        * (pathBetaO G hiso hfn S Fs hFs f) ^ ((P.degFr g : ℕ+) : ℕ) := by
  have h := rem272BetaO_comp (biratCoreOf G hfn) (fun W => birat_isOfIsotropicType hiso W)
      (fun {_ _} q => biratMet G q) hfn
      (biratFrobSection_isFrobeniusSection G hiso S hFs)
      X.path.ref_mem Y.path.ref_mem Z.path.ref_mem (pathW G hiso f) (pathW G hiso g)
      (pathGamma G hiso hfn S Fs hFs f g) (pathGamma_spec G hiso hfn S Fs hFs f g)
  refine (congrArg (rem272BetaO (biratCoreOf G hfn) (fun W => birat_isOfIsotropicType hiso W)
      (fun {_ _} q => biratMet G q) (biratFrobSection_isFrobeniusSection G hiso S hFs)
      X.path.ref_mem Z.path.ref_mem) (pathW_comp G hiso f g)).trans (h.trans ?_)
  exact congrArg (fun k : ℕ => pathGamma G hiso hfn S Fs hFs f g
      * (pathBetaO G hiso hfn S Fs hFs f) ^ k)
    (congrArg (fun q : ℕ+ => (q : ℕ)) (pathW_deg G hiso g))

/-- ★★★★★**`u` の合成則** —— model Frobenioid の合成規則そのもの。

原文 (FrdI p.102):
> the ﬁrst three entries of the data that deﬁne a morphism of C; for the ﬁnal entry, it
-/
theorem pathU_comp {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs)
    {X Y Z : PathCat S} (f : X ⟶ Y) (g : Y ⟶ Z) :
    pathU R hiso hfn S Fs hFs (f ≫ g)
      = R.bmon.map (P.Base f) (pathU R hiso hfn S Fs hFs g)
        + ((P.degFr g : ℕ+) : ℕ) • pathU R hiso hfn S Fs hFs f := by
  have hk : R.kappa X.path.ref (Additive.ofMul (pathGamma G hiso hfn S Fs hFs f g))
      = R.bmon.map (pathTheta f)
          (R.kappa Y.path.ref (Additive.ofMul (pathBetaO G hiso hfn S Fs hFs g))) :=
    R.kappa_pull (pathAf G hiso S f) (pathTheta f) (pathBetaO G hiso hfn S Fs hFs g)
      (pathGamma G hiso hfn S Fs hFs f g) (pathAf_isPullBack G hiso S f)
      (pathAf_base G hiso S f) (pathGamma_spec G hiso hfn S Fs hFs f g)
  rw [pathU_eq, pathU_eq, pathU_eq, pathBetaO_comp G hiso hfn S Fs hFs f g,
    show (Additive.ofMul (pathGamma G hiso hfn S Fs hFs f g
        * (pathBetaO G hiso hfn S Fs hFs f) ^ ((P.degFr g : ℕ+) : ℕ)))
      = Additive.ofMul (pathGamma G hiso hfn S Fs hFs f g)
        + ((P.degFr g : ℕ+) : ℕ) • Additive.ofMul (pathBetaO G hiso hfn S Fs hFs f) from rfl,
    map_add, map_nsmul, map_add, map_nsmul, hk, ← R.bmon.map_comp, pathTheta_spec,
    R.bmon.map_comp]

/-! ## ★6. 恒等射 -/

/-- ★恒等射の中央項は 1。 -/
theorem rem272BetaO_id {S : BaseSection P} (Fc : FrobenioidCore P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hmet : ∀ {A B : C} (f : A ⟶ B), IsIsometric P f)
    {Fs : ℕ+ →* SectionEnd S} (hFs : IsFrobeniusSection S Fs)
    {A : C} (hA : S.objP A) :
    rem272BetaO Fc hiso hmet hFs hA hA (𝟙 A) = 1 := by
  refine Subtype.ext ?_
  refine (rem272Beta_uniq Fc hiso hmet hFs hA hA (𝟙 A) (β := 𝟙 A)
    (n := 1) rfl ⟨P.degFr_id A, ?_⟩ (S.lift_homP hA hA (P.Base (𝟙 A))) ?_).symm
  · show IsIso (P.Base (𝟙 A))
    rw [P.Base_id]; infer_instance
  · rw [show (((Fs 1).app ⟨A, hA⟩ : End A) : A ⟶ A) = 𝟙 A from by rw [map_one]; rfl,
      S.lift_id hA, Category.id_comp, Category.id_comp]

theorem pathBetaO_id (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs)
    (X : PathCat S) : pathBetaO G hiso hfn S Fs hFs (𝟙 X) = 1 :=
  (congrArg (rem272BetaO (biratCoreOf G hfn) (fun W => birat_isOfIsotropicType hiso W)
      (fun {_ _} q => biratMet G q) (biratFrobSection_isFrobeniusSection G hiso S hFs)
      X.path.ref_mem X.path.ref_mem) (pathW_id G hiso X)).trans
    (rem272BetaO_id (biratCoreOf G hfn) (fun W => birat_isOfIsotropicType hiso W)
      (fun {_ _} q => biratMet G q) (biratFrobSection_isFrobeniusSection G hiso S hFs)
      X.path.ref_mem)

theorem pathU_id {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs)
    (X : PathCat S) : pathU R hiso hfn S Fs hFs (𝟙 X) = 0 := by
  rw [pathU_eq, pathBetaO_id]
  show R.bmon.map X.path.piBase (R.kappa X.path.ref 0) = 0
  rw [map_zero, map_zero]

/-! ## ★7. 関手 -/

/-- ★★**model Frobenioid の対象への対応** `X ↦ (Base X, cls X)`。 -/
noncomputable def pathObj {G : Frobenioid P} (R : RatFnData P G)
    {S : BaseSection P} (X : PathCat S) : ModelData.Obj R.model :=
  ⟨(P.toElem.obj X.obj).base, X.path.cls⟩

/-- ★★★★★**`𝒞̃ ⥤ 𝒞^model`** —— `Theorem 5.2, (iv)` の関手。

原文 (FrdI p.102):
> verify that these assignments determine a functor
-/
noncomputable def pathToModel {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs) :
    PathCat S ⥤ ModelData.Obj R.model where
  obj X := pathObj R X
  map {X Y} f :=
    { base := P.Base f
      div := P.Div f
      deg := P.degFr f
      u := pathU R hiso hfn S Fs hFs f
      cond := by
        show ((P.degFr f : ℕ+) : ℕ) • X.path.cls + toGpHom _ (P.Div f)
          = Φ.gpMapOn (P.Base f) Y.path.cls + R.divB _ (pathU R hiso hfn S Fs hFs f)
        rw [pathU_divB]
        abel }
  map_id X := by
    refine ModelData.Hom.ext ?_ ?_ ?_ ?_
    · exact P.Base_id X.obj
    · exact P.Div_id X.obj
    · exact P.degFr_id X.obj
    · exact pathU_id R hiso hfn S Fs hFs X
  map_comp {X Y Z} f g := by
    refine ModelData.Hom.ext ?_ ?_ ?_ ?_
    · exact P.Base_comp f g
    · exact P.Div_comp f g
    · exact P.degFr_comp f g
    · exact pathU_comp R hiso hfn S Fs hFs f g

/-- ★★★**`𝔽_Φ` への関手との 1-両立性** —— 実は**等式**で成り立つ。

★対象の底も射の 3 成分も**定義どおり一致**するので、自然同型は恒等でよい。 -/
theorem pathToModel_toElem {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs) :
    pathToModel R hiso hfn S Fs hFs ⋙ ModelData.toElem = pathForget S ⋙ P.toElem := rfl

/-! ### ★出典の紐付け(`.src`) -/

/-- ★locator —— `Theorem 5.2, (iv)` の関手 `𝒞̃ ⥤ 𝒞^model`。 -/
def pathToModel.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 102,
    item := "Theorem 5.2, (iv) — 𝒞̃ ⥤ 𝒞^model の関手",
    sectionId := "frdi-thm-5-2" }

end ABC3.Found.FrdI
