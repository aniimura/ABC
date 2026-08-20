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
  /-- `B` は **group-like**(原文: "B : D -> Mon a group-like monoid on D")。 -/
  bneg : ∀ A : D, bmon.val A → bmon.val A
  bneg_add : ∀ (A : D) (x : bmon.val A), bneg A x + x = 0

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

/-! ## ★8. 本質的全射性と忠実性 -/

/-- ★★`path` の類は、`ref` 側から見た span の類を `piBase` で運んだもの。

★これが `Theorem 5.1, (i)` の `PicObj.HasCls` と `FPPath.cls` を繋ぐ。 -/
theorem FPPath.cls_eq {S : BaseSection P} {X : C} (p : FPPath S X) :
    p.cls = Φ.gpMapOn p.piBase (spanCls p.toRef p.toRef_preStep.2 p.toObj) := by
  haveI h1 : IsIso (P.Base p.toObj) := p.toObj_preStep.2
  haveI h2 : IsIso (P.Base p.toRef) := p.toRef_preStep.2
  have e1 : p.cls = Φ.gpMapOn (inv (P.Base p.toObj))
      (toGp _ (P.Div p.toObj) - toGp _ (P.Div p.toRef)) := by
    rw [FPPath.cls, spanCls_eq]
    show -(Φ.gpMapOn (inv (P.Base p.toObj)))
      (toGp _ (P.Div p.toRef) - toGp _ (P.Div p.toObj)) = _
    rw [← map_neg, neg_sub]
  have e2 : spanCls p.toRef p.toRef_preStep.2 p.toObj
      = Φ.gpMapOn (inv (P.Base p.toRef))
        (toGp _ (P.Div p.toObj) - toGp _ (P.Div p.toRef)) := by
    rw [spanCls_eq]; rfl
  have e3 : Φ.gpMapOn (P.Base p.toRef) (Φ.gpMapOn (inv (P.Base p.toRef))
        (toGp _ (P.Div p.toObj) - toGp _ (P.Div p.toRef)))
      = toGp _ (P.Div p.toObj) - toGp _ (P.Div p.toRef) := by
    rw [← Φ.gpMapOn_comp, IsIso.hom_inv_id, Φ.gpMapOn_id]
  rw [e1, e2, FPPath.piBase, Φ.gpMapOn_comp, e3]

namespace ModelData

/-- ★底の同型から model Frobenioid の同型を作る。 -/
def isoOfBase {M : ModelData.{v, u, w} D} {A B : ModelData.Obj M}
    (θ : A.base ≅ B.base) (h : A.cls = M.phi.gpMapOn θ.hom B.cls) : A ≅ B where
  hom :=
    { base := θ.hom, div := 0, deg := 1, u := 0
      cond := by simpa using h }
  inv :=
    { base := θ.inv, div := 0, deg := 1, u := 0
      cond := by
        have : M.phi.gpMapOn θ.inv A.cls = B.cls := by
          rw [h, ← M.phi.gpMapOn_comp, θ.inv_hom_id, M.phi.gpMapOn_id]
        simpa using this.symm }
  hom_inv_id := by
    refine ModelData.Hom.ext ?_ ?_ ?_ ?_
    · exact θ.hom_inv_id
    · show M.phi.map θ.hom 0 + ((1 : ℕ+) : ℕ) • (0 : M.phi.val A.base) = 0
      simp
    · show (1 : ℕ+) * 1 = 1
      simp
    · show M.bmon.map θ.hom 0 + ((1 : ℕ+) : ℕ) • (0 : M.bmon.val A.base) = 0
      simp
  inv_hom_id := by
    refine ModelData.Hom.ext ?_ ?_ ?_ ?_
    · exact θ.inv_hom_id
    · show M.phi.map θ.inv 0 + ((1 : ℕ+) : ℕ) • (0 : M.phi.val B.base) = 0
      simp
    · show (1 : ℕ+) * 1 = 1
      simp
    · show M.bmon.map θ.inv 0 + ((1 : ℕ+) : ℕ) • (0 : M.bmon.val B.base) = 0
      simp

end ModelData

/-- ★★★★**本質的全射性** —— `Theorem 5.1, (i)` の全射性から。

原文 (FrdI p.103):
> C → C is faithful. Moreover, this functor C → C is manifestly essentially surjective [cf.
-/
theorem pathToModel_essSurj {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs) :
    (pathToModel R hiso hfn S Fs hFs).EssSurj := by
  refine ⟨fun M => ?_⟩
  obtain ⟨A, hA, ⟨e⟩⟩ := S.essSurjP M.base
  obtain ⟨Z, X₀, φ, ψ, hsφ, hsψ, hb, hspan⟩ :=
    PicObj.exists_hasCls G (A := A) (Φ.gpMapOn e.hom M.cls)
  let p : FPPath S Z.obj := ⟨A, hA, X₀, φ, ψ, hsφ, hsψ⟩
  haveI hpi : IsIso p.piBase := FPPath.isIso_piBase p
  refine ⟨⟨Z.obj, p⟩, ⟨ModelData.isoOfBase (M := R.model)
    (A := pathObj R (⟨Z.obj, p⟩ : PathCat S)) (B := M)
    ((@asIso _ _ _ _ p.piBase hpi) ≪≫ e) ?_⟩⟩
  show p.cls = Φ.gpMapOn (p.piBase ≫ e.hom) M.cls
  rw [Φ.gpMapOn_comp, ← hspan]
  exact p.cls_eq

/-- ★同型に沿った `MonoidOn.map` は単射。 -/
theorem MonoidOn.map_injective_of_isIso (M : MonoidOn.{v, u, w} D) {A B : D} (θ : A ⟶ B)
    [IsIso θ] : Function.Injective (M.map θ) := by
  intro x y hxy
  have h : M.map (inv θ) (M.map θ x) = M.map (inv θ) (M.map θ y) := congrArg _ hxy
  rwa [← M.map_comp, IsIso.inv_hom_id, M.map_id, ← M.map_comp, IsIso.inv_hom_id,
    M.map_id] at h

/-- ★★★**3 分解の 3 成分が一致すれば射も一致する**(`Remark 2.7.2` の一意性の裏返し)。 -/
theorem eq_of_rem272Beta_eq {S : BaseSection P} (Fc : FrobenioidCore P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hmet : ∀ {A B : C} (f : A ⟶ B), IsIsometric P f)
    {Fs : ℕ+ →* SectionEnd S} (hFs : IsFrobeniusSection S Fs)
    {A B : C} (hA : S.objP A) (hB : S.objP B) (φ φ' : A ⟶ B)
    (hd : P.degFr φ = P.degFr φ') (hb : P.Base φ = P.Base φ')
    (hbeta : rem272Beta Fc hiso hmet hFs hA hB φ = rem272Beta Fc hiso hmet hFs hA hB φ') :
    φ = φ' := by
  conv_lhs => rw [rem272Beta_spec Fc hiso hmet hFs hA hB φ]
  rw [hd, hb, hbeta, ← rem272Beta_spec Fc hiso hmet hFs hA hB φ']

/-- ★★★★**忠実性** —— `u_f` から `β_f` が復元され、3 分解の一意性で `w_f` が定まる。

原文 (FrdI p.103):
> C → C is faithful. Moreover, this functor C → C is manifestly essentially surjective [cf.
-/
theorem pathToModel_faithful {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs) :
    (pathToModel R hiso hfn S Fs hFs).Faithful := by
  refine ⟨fun {X Y} {f g} h => ?_⟩
  have hbase : P.Base f = P.Base g := congrArg ModelData.Hom.base h
  have hdegf : P.degFr f = P.degFr g := congrArg ModelData.Hom.deg h
  have hu : pathU R hiso hfn S Fs hFs f = pathU R hiso hfn S Fs hFs g :=
    congrArg ModelData.Hom.u h
  haveI hpi : IsIso X.path.piBase := FPPath.isIso_piBase X.path
  have hbeta : pathBetaO G hiso hfn S Fs hFs f = pathBetaO G hiso hfn S Fs hFs g :=
    Additive.ofMul.injective ((R.kappa X.path.ref).injective
      (R.bmon.map_injective_of_isIso X.path.piBase hu))
  have hth : pathTheta f = pathTheta g := by rw [pathTheta, pathTheta, hbase]
  have hw : pathW G hiso f = pathW G hiso g :=
    eq_of_rem272Beta_eq (biratCoreOf G hfn) (fun W => birat_isOfIsotropicType hiso W)
      (fun {_ _} q => biratMet G q) (biratFrobSection_isFrobeniusSection G hiso S hFs)
      X.path.ref_mem Y.path.ref_mem (pathW G hiso f) (pathW G hiso g)
      (by show biratDeg (pathW G hiso f) = biratDeg (pathW G hiso g)
          rw [pathW_deg, pathW_deg, hdegf])
      (by show biratBase (pathW G hiso f) = biratBase (pathW G hiso g)
          rw [biratBase_pathW, biratBase_pathW, hth])
      (congrArg Subtype.val hbeta)
  haveI hpiY : IsIso (Y.path.piHom G hiso) := Y.path.isIso_piHom G hiso
  have hfb : (toBiratCat P G).map f = (toBiratCat P G).map g := by
    have e1 := pathW_spec G hiso f
    have e2 := pathW_spec G hiso g
    rw [hw] at e1
    exact (cancel_mono (Y.path.piHom G hiso)).mp (e1.symm.trans e2)
  haveI := toBiratCat_faithful (P := P) (G := G)
  exact (toBiratCat P G).map_injective hfb

/-! ## ★9. 充満性

原文 (FrdI p.103):
> C → C is faithful. Moreover, this functor C → C is manifestly essentially surjective [cf.

★★充満性は 2 段でできている:
1. `Theorem 5.1, (ii)` で **(Base, Div, degFr) の合う射 `f₀`** を取る。
2. `u` のずれ `v ∈ ker(Div_B)` を、`Definition 1.3, (iv)(a)` の 3 分解
   `f₀ = γ ≫ β₀ ≫ α` の **`γ` の終域に単元を挿す**ことで吸収する。

★★2 が要点である。単元を**左**から掛けると `u ↦ u + deg(f)·u_ε`、**右**から掛けると
`u ↦ u + B(Base f)(u_ε)` で、どちらも一般には `ker(Div_B)` 全体を走らない。
★しかし `γ` は Frobenius 型なので **`Base γ` は同型**であり、`B(Base γ)` は全単射。
したがって中間に挿せば `ker(Div_B)` の任意の元が実現できる。 -/

/-- ★`𝒞` の単元の像は `𝒞^birat` の単元。 -/
theorem otimes_toBirat (G : Frobenioid P) {A : C} (ε : OTimes P A) :
    ((toBiratCat P G).map ((ε : End A) : A ⟶ A) : End ((toBiratCat P G).obj A))
      ∈ OTimes (biratPre P G) ((toBiratCat P G).obj A) := by
  haveI : IsIso ((ε : End A) : A ⟶ A) := (isUnit_iff_isIso (ε : End A)).mp ε.2.2
  refine ⟨⟨?_, ?_⟩, ?_⟩
  · show biratBase ((toBiratCat P G).map ((ε : End A) : A ⟶ A))
      = biratBase (𝟙 ((toBiratCat P G).obj A))
    rw [biratBase_toBiratMap,
      show (𝟙 ((toBiratCat P G).obj A)) = (toBiratCat P G).map (𝟙 A) from
        ((toBiratCat P G).map_id A).symm, biratBase_toBiratMap]
    exact ε.2.1.1
  · show biratDeg ((toBiratCat P G).map ((ε : End A) : A ⟶ A)) = 1
    rw [biratDeg_toBiratMap]
    exact ε.2.1.2
  · refine (isUnit_iff_isIso (((toBiratCat P G).map ((ε : End A) : A ⟶ A)) :
      End ((toBiratCat P G).obj A))).mpr ?_
    exact (toBiratCat P G).map_isIso _

/-- ★単元の `w` は底が恒等。 -/
theorem pathW_otimes_base (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {S : BaseSection P} {X : PathCat S} (ε : OTimes P X.obj) :
    biratBase (pathW G hiso ((ε : End X.obj) : X.obj ⟶ X.obj)) = 𝟙 _ := by
  haveI : IsIso X.path.piBase := FPPath.isIso_piBase X.path
  rw [biratBase_pathW, pathTheta,
    show P.Base ((ε : End X.obj) : X.obj ⟶ X.obj) = 𝟙 _ from by
      have h : P.Base ((ε : End X.obj) : X.obj ⟶ X.obj) = P.Base (𝟙 X.obj) := ε.2.1.1
      rwa [P.Base_id] at h,
    Category.id_comp, IsIso.inv_hom_id]
  rfl

theorem pathW_otimes_deg (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {S : BaseSection P} {X : PathCat S} (ε : OTimes P X.obj) :
    biratDeg (pathW G hiso ((ε : End X.obj) : X.obj ⟶ X.obj)) = 1 := by
  rw [pathW_deg]
  exact ε.2.1.2

/-- ★★単元の 3 分解の中央項は `w` そのもの。 -/
theorem pathBetaO_otimes (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs)
    {X : PathCat S} (ε : OTimes P X.obj) :
    ((pathBetaO G hiso hfn S Fs hFs ((ε : End X.obj) : X.obj ⟶ X.obj) :
        End ((toBiratCat P G).obj X.path.ref)) : _ ⟶ _)
      = pathW G hiso ((ε : End X.obj) : X.obj ⟶ X.obj) := by
  set S' := biratBaseSection G hiso S with hS'
  set w := pathW G hiso ((ε : End X.obj) : X.obj ⟶ X.obj) with hw
  have hbw : biratBase w = 𝟙 _ := pathW_otimes_base G hiso ε
  have hdw : biratDeg w = 1 := pathW_otimes_deg G hiso ε
  have hbid : IsBaseIdentity (biratPre P G) w := by
    show biratBase w = (biratPre P G).Base (𝟙 _)
    rw [hbw, (biratPre P G).Base_id]
    rfl
  have hps : IsPreStep (biratPre P G) w := ⟨hdw, by
    show IsIso (biratBase w); rw [hbw]; infer_instance⟩
  have hlift : S'.lift X.path.ref_mem X.path.ref_mem ((biratPre P G).Base w)
      = 𝟙 ((toBiratCat P G).obj X.path.ref) :=
    (congrArg (S'.lift X.path.ref_mem X.path.ref_mem) hbid).trans (S'.lift_id X.path.ref_mem)
  have hF1 : ((((biratFrobSection G hiso S Fs) 1).app ⟨(toBiratCat P G).obj X.path.ref,
      X.path.ref_mem⟩ : End _) : _ ⟶ _) = 𝟙 ((toBiratCat P G).obj X.path.ref) := by
    rw [map_one]; rfl
  have htriv : (𝟙 ((toBiratCat P G).obj X.path.ref)) ≫ w
      ≫ (𝟙 ((toBiratCat P G).obj X.path.ref)) = w := by
    rw [Category.id_comp, Category.comp_id]
  have step1 : (𝟙 ((toBiratCat P G).obj X.path.ref)) ≫ w
        ≫ (𝟙 ((toBiratCat P G).obj X.path.ref))
      = (𝟙 ((toBiratCat P G).obj X.path.ref)) ≫ w
        ≫ S'.lift X.path.ref_mem X.path.ref_mem ((biratPre P G).Base w) :=
    congrArg (fun t : (toBiratCat P G).obj X.path.ref ⟶ (toBiratCat P G).obj X.path.ref =>
      (𝟙 ((toBiratCat P G).obj X.path.ref)) ≫ w ≫ t) hlift.symm
  have step2 : (𝟙 ((toBiratCat P G).obj X.path.ref)) ≫ w
        ≫ S'.lift X.path.ref_mem X.path.ref_mem ((biratPre P G).Base w)
      = ((((biratFrobSection G hiso S Fs) 1).app ⟨(toBiratCat P G).obj X.path.ref,
          X.path.ref_mem⟩ : End _) : _ ⟶ _) ≫ w
        ≫ S'.lift X.path.ref_mem X.path.ref_mem ((biratPre P G).Base w) :=
    congrArg (fun t : (toBiratCat P G).obj X.path.ref ⟶ (toBiratCat P G).obj X.path.ref =>
      t ≫ w ≫ S'.lift X.path.ref_mem X.path.ref_mem ((biratPre P G).Base w)) hF1.symm
  exact (rem272Beta_uniq (biratCoreOf G hfn) (fun W => birat_isOfIsotropicType hiso W)
    (fun {_ _} q => biratMet G q) (biratFrobSection_isFrobeniusSection G hiso S hFs)
    X.path.ref_mem X.path.ref_mem w (n := 1) hbid hps
    (S'.lift_homP X.path.ref_mem X.path.ref_mem ((biratPre P G).Base w))
    (htriv.symm.trans (step1.trans step2))).symm

/-- ★★★**単元の `u` は `κ` そのもの** —— `path` の同型に沿った `κ` の自然性から。 -/
theorem pathU_otimes {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs)
    {X : PathCat S} (ε : OTimes P X.obj) :
    pathU R hiso hfn S Fs hFs ((ε : End X.obj) : X.obj ⟶ X.obj)
      = R.kappa X.obj (Additive.ofMul
          ⟨(toBiratCat P G).map ((ε : End X.obj) : X.obj ⟶ X.obj), otimes_toBirat G ε⟩) := by
  haveI hpi : IsIso (X.path.piHom G hiso) := X.path.isIso_piHom G hiso
  have hsq : X.path.piHom G hiso
        ≫ ((pathBetaO G hiso hfn S Fs hFs ((ε : End X.obj) : X.obj ⟶ X.obj) :
            End ((toBiratCat P G).obj X.path.ref)) : _ ⟶ _)
      = ((⟨(toBiratCat P G).map ((ε : End X.obj) : X.obj ⟶ X.obj), otimes_toBirat G ε⟩ :
            OTimes (biratPre P G) ((toBiratCat P G).obj X.obj)) :
            End ((toBiratCat P G).obj X.obj)) ≫ X.path.piHom G hiso := by
    rw [pathBetaO_otimes G hiso hfn S Fs hFs ε]
    exact pathW_spec G hiso ((ε : End X.obj) : X.obj ⟶ X.obj)
  refine (pathU_eq R hiso hfn S Fs hFs _).trans ?_
  exact (R.kappa_pull (X.path.piHom G hiso) X.path.piBase
    (pathBetaO G hiso hfn S Fs hFs ((ε : End X.obj) : X.obj ⟶ X.obj))
    ⟨(toBiratCat P G).map ((ε : End X.obj) : X.obj ⟶ X.obj), otimes_toBirat G ε⟩
    (isPullBack_of_isIso (biratPre P G) (X.path.piHom G hiso))
    (FPPath.biratBase_piHom G hiso X.path) hsq).symm

/-! ### ★`Theorem 5.1, (ii)` を当てるための道具 -/

/-- ★★`path` が `ref` の側で定める類(`Theorem 5.1, (i)` の `PicObj.HasCls` の類)。 -/
noncomputable def FPPath.refCls {S : BaseSection P} {X : C} (p : FPPath S X) :
    Gp (Φ.val (P.toElem.obj p.ref).base) :=
  spanCls p.toRef p.toRef_preStep.2 p.toObj

/-- ★path が定める `PicObj`。 -/
noncomputable def FPPath.picObj {S : BaseSection P} {X : C} (p : FPPath S X) :
    PicObj P p.ref :=
  ⟨X, @asIso _ _ _ _ p.piBase (FPPath.isIso_piBase p)⟩

theorem FPPath.picObj_hasCls {S : BaseSection P} {X : C} (p : FPPath S X) :
    p.picObj.HasCls (spanCls p.toRef p.toRef_preStep.2 p.toObj) := by
  haveI h1 : IsIso (P.Base p.toObj) := p.toObj_preStep.2
  refine ⟨p.vertex, p.toRef, p.toObj, p.toRef_preStep, p.toObj_preStep, ?_, rfl⟩
  show P.Base p.toObj ≫ p.piBase = P.Base p.toRef
  rw [FPPath.piBase, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp]

/-- ★★`path` が `ref` の側で定める類。 -/
theorem FPPath.cls_eq' {S : BaseSection P} {X : C} (p : FPPath S X) :
    p.cls = Φ.gpMapOn p.piBase p.refCls := p.cls_eq

theorem FPPath.picObj_hasCls' {S : BaseSection P} {X : C} (p : FPPath S X) :
    p.picObj.HasCls p.refCls := p.picObj_hasCls

/-- ★★★★**model の条件式 ⟹ `Theorem 5.1, (ii)` の所属条件の中身**。 -/
theorem model_cond_eq {G : Frobenioid P} (R : RatFnData P G)
    {S : BaseSection P} {X Y : PathCat S} (n : ℕ+)
    (θ₀ : (P.toElem.obj X.obj).base ⟶ (P.toElem.obj Y.obj).base)
    (dv : Φ.val (P.toElem.obj X.obj).base) (uu : R.bmon.val (P.toElem.obj X.obj).base)
    (hcond : ((n : ℕ+) : ℕ) • X.path.cls + toGp _ dv
      = Φ.gpMapOn θ₀ Y.path.cls + R.divB _ uu) :
    ((n : ℕ+) : ℕ) • X.path.refCls + toGp _ (Φ.map (inv X.path.piBase) dv)
        - Φ.gpMapOn ((inv X.path.piBase ≫ θ₀) ≫ Y.path.piBase) Y.path.refCls
      = R.divB _ (R.bmon.map (inv X.path.piBase) uu) := by
  haveI hiX : IsIso X.path.piBase := FPPath.isIso_piBase X.path
  haveI hiY : IsIso Y.path.piBase := FPPath.isIso_piBase Y.path
  have hcX : Φ.gpMapOn (inv X.path.piBase) X.path.cls = X.path.refCls := by
    rw [X.path.cls_eq', ← Φ.gpMapOn_comp, IsIso.inv_hom_id, Φ.gpMapOn_id]
  have hcY : Φ.gpMapOn (inv X.path.piBase) (Φ.gpMapOn θ₀ Y.path.cls)
      = Φ.gpMapOn ((inv X.path.piBase ≫ θ₀) ≫ Y.path.piBase) Y.path.refCls := by
    rw [Y.path.cls_eq', ← Φ.gpMapOn_comp, ← Φ.gpMapOn_comp]
  have h := congrArg (⇑(Φ.gpMapOn (inv X.path.piBase))) hcond
  rw [map_add, map_add, map_nsmul, hcX, hcY, ← R.divB_nat,
    show Φ.gpMapOn (inv X.path.piBase) (toGp (Φ.val (P.toElem.obj X.obj).base) dv)
      = toGp _ (Φ.map (inv X.path.piBase) dv) from gpMap_toGp _ _ _] at h
  rw [h]
  abel

/-- ★`Div_B` の像は `Φ^birat` に入る(`κ` が全単射だから)。 -/
theorem divB_mem_phiBiratAt {G : Frobenioid P} (R : RatFnData P G) {A : C}
    (y : R.bmon.val (P.toElem.obj A).base) :
    R.divB _ y ∈ phiBiratAt P G ((toBiratCat P G).obj A) := by
  refine ⟨((Additive.toMul ((R.kappa A).symm y) :
      OTimes (biratPre P G) ((toBiratCat P G).obj A)) : End _),
    (Additive.toMul ((R.kappa A).symm y)).2, ?_⟩
  refine (R.kappa_divB A (Additive.toMul ((R.kappa A).symm y))).symm.trans ?_
  exact congrArg _ ((R.kappa A).apply_symm_apply y)

/-- ★★★**`Theorem 5.1, (ii)` から (Base, Div, degFr) の合う射を取る**。 -/
theorem exists_hom_of_model {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y) {S : BaseSection P} {X Y : PathCat S} (n : ℕ+)
    (θ₀ : (P.toElem.obj X.obj).base ⟶ (P.toElem.obj Y.obj).base)
    (dv : Φ.val (P.toElem.obj X.obj).base) (uu : R.bmon.val (P.toElem.obj X.obj).base)
    (hcond : ((n : ℕ+) : ℕ) • X.path.cls + toGp _ dv
      = Φ.gpMapOn θ₀ Y.path.cls + R.divB _ uu) :
    ∃ f₀ : X.obj ⟶ Y.obj, P.degFr f₀ = n ∧ P.Base f₀ = θ₀ ∧ P.Div f₀ = dv := by
  haveI hiX : IsIso X.path.piBase := FPPath.isIso_piBase X.path
  haveI hiY : IsIso Y.path.piBase := FPPath.isIso_piBase Y.path
  have hmem := divB_mem_phiBiratAt R (R.bmon.map (inv X.path.piBase) uu)
  rw [← model_cond_eq R n θ₀ dv uu hcond] at hmem
  obtain ⟨f₀, hdeg0, hbase0, hdiv0⟩ :=
    (thm_5_1_ii G hiso (S.frobTrivial X.path.ref_mem) (S.frobTrivial Y.path.ref_mem)
      X.path.picObj Y.path.picObj X.path.picObj_hasCls' Y.path.picObj_hasCls'
      ((inv X.path.piBase ≫ θ₀) ≫ Y.path.piBase) n dv).mpr hmem
  refine ⟨f₀, hdeg0, hbase0.trans ?_, hdiv0⟩
  show X.path.piBase ≫ ((inv X.path.piBase ≫ θ₀) ≫ Y.path.piBase) ≫ inv Y.path.piBase = θ₀
  rw [Category.assoc, IsIso.hom_inv_id, Category.comp_id, ← Category.assoc,
    IsIso.hom_inv_id, Category.id_comp]

/-! ### ★`𝒞` の合成で書いた `pathU` -/

/-- ★`pathU` を **`𝒞` の射**で書いた版(`PathCat` の合成と構文が違うため)。 -/
noncomputable def pathUC {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs)
    (X Y : PathCat S) (f : X.obj ⟶ Y.obj) : R.bmon.val (P.toElem.obj X.obj).base :=
  pathU R hiso hfn S Fs hFs (X := X) (Y := Y) f

theorem pathUC_comp {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs)
    (X Y Z : PathCat S) (f : X.obj ⟶ Y.obj) (g : Y.obj ⟶ Z.obj) :
    pathUC R hiso hfn S Fs hFs X Z (f ≫ g)
      = R.bmon.map (P.Base f) (pathUC R hiso hfn S Fs hFs Y Z g)
        + ((P.degFr g : ℕ+) : ℕ) • pathUC R hiso hfn S Fs hFs X Y f :=
  pathU_comp R hiso hfn S Fs hFs (X := X) (Y := Y) (Z := Z) f g

theorem pathUC_otimes {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs)
    (X : PathCat S) (ε : OTimes P X.obj) :
    pathUC R hiso hfn S Fs hFs X X ((ε : End X.obj) : X.obj ⟶ X.obj)
      = R.kappa X.obj (Additive.ofMul
          ⟨(toBiratCat P G).map ((ε : End X.obj) : X.obj ⟶ X.obj), otimes_toBirat G ε⟩) :=
  pathU_otimes R hiso hfn S Fs hFs ε

/-- ★★★**単元を中間に挿すと `u` がちょうど 1 回ぶんずれる**。 -/
theorem pathUC_insert {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs)
    (X XU Y : PathCat S) (γ : X.obj ⟶ XU.obj) (h : XU.obj ⟶ Y.obj) (hdh : P.degFr h = 1)
    (δ : OTimes P XU.obj) :
    pathUC R hiso hfn S Fs hFs X Y (γ ≫ ((δ : End XU.obj) : XU.obj ⟶ XU.obj) ≫ h)
      = pathUC R hiso hfn S Fs hFs X Y (γ ≫ h)
        + R.bmon.map (P.Base γ) (R.kappa XU.obj (Additive.ofMul
            ⟨(toBiratCat P G).map ((δ : End XU.obj) : XU.obj ⟶ XU.obj),
              otimes_toBirat G δ⟩)) := by
  have hbδ : P.Base ((δ : End XU.obj) : XU.obj ⟶ XU.obj) = 𝟙 _ := by
    have hh : P.Base ((δ : End XU.obj) : XU.obj ⟶ XU.obj) = P.Base (𝟙 XU.obj) := δ.2.1.1
    rwa [P.Base_id] at hh
  have hdδ : P.degFr ((δ : End XU.obj) : XU.obj ⟶ XU.obj) = 1 := δ.2.1.2
  have hdc : P.degFr (((δ : End XU.obj) : XU.obj ⟶ XU.obj) ≫ h) = 1 := by
    rw [P.degFr_comp, hdh, hdδ, mul_one]
  have e1 : pathUC R hiso hfn S Fs hFs XU Y (((δ : End XU.obj) : XU.obj ⟶ XU.obj) ≫ h)
      = pathUC R hiso hfn S Fs hFs XU Y h
        + R.kappa XU.obj (Additive.ofMul
            ⟨(toBiratCat P G).map ((δ : End XU.obj) : XU.obj ⟶ XU.obj),
              otimes_toBirat G δ⟩) := by
    rw [pathUC_comp R hiso hfn S Fs hFs XU XU Y ((δ : End XU.obj) : XU.obj ⟶ XU.obj) h,
      hbδ, hdh, pathUC_otimes R hiso hfn S Fs hFs XU δ,
      show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul]
    exact congrArg (fun t => t + _) (R.bmon.map_id _ _)
  rw [pathUC_comp R hiso hfn S Fs hFs X XU Y γ
      (((δ : End XU.obj) : XU.obj ⟶ XU.obj) ≫ h),
    pathUC_comp R hiso hfn S Fs hFs X XU Y γ h,
    hdc, hdh, e1, map_add, show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul]
  abel

/-- ★★★★★**`u` を `ker(Div_B)` のぶんだけ自由にずらせる**。

★`Definition 1.3, (iv)(a)` の 3 分解 `f₀ = γ ≫ β₀ ≫ α` の
**Frobenius 部分 `γ` の終域**に単元を挿す。★そこの底 `Base γ` は**同型**なので
`B(Base γ)` は全単射で、`ker(Div_B)` の任意の元が実現できる。 -/
theorem exists_hom_shift {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs)
    (X Y : PathCat S) (f₀ : X.obj ⟶ Y.obj)
    (v : R.bmon.val (P.toElem.obj X.obj).base) (hv : R.divB _ v = 0) :
    ∃ f : X.obj ⟶ Y.obj, P.degFr f = P.degFr f₀ ∧ P.Div f = P.Div f₀
      ∧ P.Base f = P.Base f₀
      ∧ pathUC R hiso hfn S Fs hFs X Y f = pathUC R hiso hfn S Fs hFs X Y f₀ + v := by
  obtain ⟨U, V, γ, β₀, αp, hfac, hγ, hβs, hαpb⟩ := G.core.arbFactor f₀
  haveI hbγ : IsIso (P.Base γ) := hγ.2
  obtain ⟨pU⟩ := exists_fpPath G S U
  have hx : R.divB _ (R.bmon.map (inv (P.Base γ)) v) = 0 := by
    rw [R.divB_nat, hv, map_zero]
  set zh : OTimes (biratPre P G) ((toBiratCat P G).obj U) :=
    Additive.toMul ((R.kappa U).symm (R.bmon.map (inv (P.Base γ)) v)) with hzh
  have hkz : R.kappa U (Additive.ofMul zh) = R.bmon.map (inv (P.Base γ)) v :=
    (R.kappa U).apply_symm_apply _
  have hzd : biratDivGp ((zh : End ((toBiratCat P G).obj U)) :
      (toBiratCat P G).obj U ⟶ (toBiratCat P G).obj U) = 0 :=
    (R.kappa_divB U zh).symm.trans ((congrArg _ hkz).trans hx)
  obtain ⟨δ, hδmem, hδeq⟩ := (phiBiratAt_ker ((toBiratCat P G).obj U) zh.2).mp hzd
  have hdh : P.degFr (β₀ ≫ αp) = 1 := by
    rw [P.degFr_comp, (G.core.pullBackLB αp hαpb).2, hβs.1, mul_one]
  have hue : IsUnitEquivalent P f₀ (γ ≫ ((δ : U ⟶ U)) ≫ (β₀ ≫ αp)) :=
    ⟨U, γ, β₀ ≫ αp, δ, hδmem, hfac, rfl⟩
  obtain ⟨hd, hz, hb⟩ := prop_3_3_ii_necessity P hue
  refine ⟨γ ≫ ((δ : U ⟶ U)) ≫ (β₀ ≫ αp), hd.symm, hz.symm, hb.symm, ?_⟩
  refine (pathUC_insert R hiso hfn S Fs hFs X (⟨U, pU⟩ : PathCat S) Y γ (β₀ ≫ αp) hdh
    ⟨(δ : U ⟶ U), hδmem⟩).trans ?_
  refine congrArg₂ (· + ·)
    (congrArg (pathUC R hiso hfn S Fs hFs X Y) hfac.symm) ?_
  have hk : R.kappa U (Additive.ofMul
      (⟨(toBiratCat P G).map ((δ : U ⟶ U)), otimes_toBirat G ⟨(δ : U ⟶ U), hδmem⟩⟩ :
        OTimes (biratPre P G) ((toBiratCat P G).obj U)))
      = R.bmon.map (inv (P.Base γ)) v :=
    (congrArg (fun t : OTimes (biratPre P G) ((toBiratCat P G).obj U) =>
      R.kappa U (Additive.ofMul t)) (Subtype.ext hδeq.symm :
        (⟨(toBiratCat P G).map ((δ : U ⟶ U)),
          otimes_toBirat G ⟨(δ : U ⟶ U), hδmem⟩⟩ :
        OTimes (biratPre P G) ((toBiratCat P G).obj U)) = zh)).trans hkz
  rw [hk, ← R.bmon.map_comp, IsIso.hom_inv_id, R.bmon.map_id]

/-- ★★★★★**model の射のデータから `𝒞` の射を作る**(充満性の中身)。 -/
theorem exists_map_eq {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs)
    (X Y : PathCat S)
    (mb : (P.toElem.obj X.obj).base ⟶ (P.toElem.obj Y.obj).base)
    (mdv : Φ.val (P.toElem.obj X.obj).base) (mdg : ℕ+)
    (mu : R.bmon.val (P.toElem.obj X.obj).base)
    (hcond : ((mdg : ℕ+) : ℕ) • X.path.cls + toGp _ mdv
      = Φ.gpMapOn mb Y.path.cls + R.divB _ mu) :
    ∃ f : X.obj ⟶ Y.obj, P.Base f = mb ∧ P.Div f = mdv ∧ P.degFr f = mdg
      ∧ pathUC R hiso hfn S Fs hFs X Y f = mu := by
  obtain ⟨f₀, hd0, hb0, hz0⟩ :=
    exists_hom_of_model R hiso (X := X) (Y := Y) mdg mb mdv mu hcond
  have hdivf : R.divB _ (pathUC R hiso hfn S Fs hFs X Y f₀)
      = R.divB ((P.toElem.obj X.obj).base) mu := by
    show R.divB _ (pathU R hiso hfn S Fs hFs (X := X) (Y := Y) f₀) = _
    rw [pathU_divB R hiso hfn S Fs hFs (X := X) (Y := Y) f₀, hd0, hb0, hz0, hcond]
    abel
  have hnegd : R.divB _ (R.bneg _ (pathUC R hiso hfn S Fs hFs X Y f₀))
      + R.divB _ (pathUC R hiso hfn S Fs hFs X Y f₀) = 0 := by
    rw [← map_add, R.bneg_add, map_zero]
  have hv : R.divB _ (R.bneg _ (pathUC R hiso hfn S Fs hFs X Y f₀) + mu) = 0 := by
    rw [map_add, ← hdivf]
    exact hnegd
  obtain ⟨f, hd, hz, hb, hu⟩ := exists_hom_shift R hiso hfn S Fs hFs X Y f₀
    (R.bneg _ (pathUC R hiso hfn S Fs hFs X Y f₀) + mu) hv
  have hcancel : pathUC R hiso hfn S Fs hFs X Y f₀
      + R.bneg _ (pathUC R hiso hfn S Fs hFs X Y f₀) = 0 := by
    rw [add_comm]; exact R.bneg_add _ _
  refine ⟨f, hb.trans hb0, hz.trans hz0, hd.trans hd0, ?_⟩
  rw [hu, ← add_assoc, hcancel, zero_add]

/-- ★★★★★**充満性**。 -/
theorem pathToModel_full {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs) :
    (pathToModel R hiso hfn S Fs hFs).Full := by
  refine ⟨fun {X Y} m => ?_⟩
  obtain ⟨f, hb, hz, hd, hu⟩ :=
    exists_map_eq R hiso hfn S Fs hFs X Y m.base m.div m.deg m.u m.cond
  exact ⟨f, ModelData.Hom.ext hb hz hd hu⟩

/-! ## ★10. `Theorem 5.2, (iv)` -/

/-- ★★★★★★**`𝒞̃ ⥤ 𝒞^model` は圏同値**。 -/
theorem pathToModel_isEquivalence {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs) :
    (pathToModel R hiso hfn S Fs hFs).IsEquivalence := by
  haveI := pathToModel_faithful R hiso hfn S Fs hFs
  haveI := pathToModel_full R hiso hfn S Fs hFs
  haveI := pathToModel_essSurj R hiso hfn S Fs hFs
  exact Functor.IsEquivalence.mk

/-- ★★★★★★★**[FrdI] Theorem 5.2, (iv)** —— `𝒞` は model Frobenioid と圏同値。

原文 (FrdI p.101):
> that is 1-compatible with the functors C → FΦ, C → FΦ.
-/
noncomputable def thm_5_2_iv {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs) :
    C ≌ ModelData.Obj R.model := by
  haveI := pathForget_isEquivalence G S
  haveI := pathToModel_isEquivalence R hiso hfn S Fs hFs
  exact (pathForget S).asEquivalence.symm.trans
    (pathToModel R hiso hfn S Fs hFs).asEquivalence

/-! ### ★出典の紐付け(`.src`) -/

/-- ★locator —— `Theorem 5.2, (iv)` の関手 `𝒞̃ ⥤ 𝒞^model`。 -/
def pathToModel.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 102,
    item := "Theorem 5.2, (iv) — 𝒞̃ ⥤ 𝒞^model の関手",
    sectionId := "frdi-thm-5-2" }

/-- ★★★locator —— **[FrdI] Theorem 5.2, (iv)** そのもの。 -/
def thm_5_2_iv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 101,
    item := "Theorem 5.2, (iv) — model 型の Frobenioid は model Frobenioid と圏同値",
    sectionId := "frdi-thm-5-2" }

end ABC3.Found.FrdI
