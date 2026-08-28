/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.TensorSurj
import ABC3.Found.Arakelov.PullbackGen
import ABC3.Meta.Claim

/-!
# 底変換の道具立て —— `𝓞_K ⊗_{𝓞_F} Γ(L) → Γ(f^*L)` を作る（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

## ★★★★★★★★これは何か

底変換 `deg_K(L|_{Spec 𝓞_K}) = deg_F(L)` の**有限側**を出すために、
切断の引き戻し `pullSec` から係数拡大の写像

    `μ_L : 𝓞_K ⊗_{𝓞_F} Γ(L) → Γ(f^*L)`

を作る。★これが同型であることが最後の目標であり、道は `§9-785` に書いた
「テンソルで回す」——`μ_A ⊗ μ_{A⁻¹} ≅ μ_{Ō}` は同型なので `μ_A` は全射、
したがって `Module.Invertible.bijective_of_surjective` で同型。

## ★★★★★本ファイルが用意するもの

| もの | 内容 |
|---|---|
| `gammaModPre` | ★`Γ(L)` を `R`-加群として見る（`Γ-Spec` 同型で係数制限） |
| `pullSecTop` | ★★`pullSec` の `⊤` 成分（型を `op ⊤` で宣言し直したもの） |
| `pullSecTop_add` / `pullSecTop_smul` | ★★★加法性と `ρ`-半線形性 |
| `isScalarTower_compHom` | ★係数制限は塔をなす |
| `bcLift` / `bcLift_tmul` | ★★★★★半線形写像から `T ⊗_S M →ₗ[T] N` を作る |

★★`pullSecTop_smul` の中身は在庫の `pullSec_smul`（`§9-741`）と
mathlib の `Scheme.ΓSpecIso_inv_naturality` の合成である。

## ★型の摩擦について（記録）

`pullSec f L ⊤` の行き先は `((f^*L).obj (op ((Opens.map f.base).obj ⊤)))` であり、
`(Opens.map f.base).obj ⊤ = ⊤` は `rfl` だが**インスタンス探索はそれを見ない**。
★そこで `pullSecTop` として**戻り値の型を `op ⊤` で宣言し直す**
——`gammaSheafifyM`（`§9-780`）と同じ手口である。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace
open ABC3.Found.GenEll
open scoped TensorProduct

/-! ## ★★★★`Γ(L)` を係数環の側で見る -/

/-- ★**前層の大域切断を `R`-加群として見る**（`Γ-Spec` 同型で係数制限）。

★`moduleSpecΓFunctor` の前層版である（`gammaEqRestrict` と同じ形）。 -/
noncomputable def gammaModPre (R : CommRingCat.{0}) (L : (Spec R).PresheafOfModules) :
    ModuleCat R :=
  (ModuleCat.restrictScalars (Scheme.ΓSpecIso R).inv.hom).obj (L.obj (op ⊤))

/-! ## ★★★★★`pullSec` の `⊤` 成分 -/

/-- ★★**切断の引き戻しの `⊤` 成分**（戻り値の型を `op ⊤` で宣言し直したもの）。 -/
noncomputable def pullSecTop {R T : CommRingCat.{0}} (h : R ⟶ T)
    (L : (Spec R).PresheafOfModules) (s : (L.obj (op ⊤) : Type)) :
    (((pullbackPre (Spec.map h)).obj L).obj (op (⊤ : (Spec T).Opens)) : Type) :=
  pullSec (Spec.map h) L ⊤ s

/-- ★加法性。 -/
theorem pullSecTop_add {R T : CommRingCat.{0}} (h : R ⟶ T)
    (L : (Spec R).PresheafOfModules) (s t : (L.obj (op ⊤) : Type)) :
    pullSecTop h L (s + t) = pullSecTop h L s + pullSecTop h L t := by
  show pullSec (Spec.map h) L ⊤ (s + t) = _
  exact map_add (((PresheafOfModules.pullbackPushforwardAdjunction
    (pullbackPhi (Spec.map h))).unit.app L).app (op ⊤)).hom s t

/-- ★★★**係数環の側で見た半線形性**。

★在庫の `pullSec_smul`（`§9-741`）と mathlib の `Scheme.ΓSpecIso_inv_naturality` の合成。 -/
theorem pullSecTop_smul {R T : CommRingCat.{0}} (h : R ⟶ T)
    (L : (Spec R).PresheafOfModules) (a : (R : Type)) (s : (L.obj (op ⊤) : Type)) :
    pullSecTop h L (((Scheme.ΓSpecIso R).inv.hom a) • s)
      = ((Scheme.ΓSpecIso T).inv.hom (h.hom a)) • pullSecTop h L s := by
  show pullSec (Spec.map h) L ⊤ (((Scheme.ΓSpecIso R).inv.hom a) • s) = _
  rw [pullSec_smul]
  congr 1
  exact (congrArg (fun (m : _ ⟶ _) => CommRingCat.Hom.hom m a)
    (Scheme.ΓSpecIso_inv_naturality h)).symm

/-! ## ★★★★★★半線形写像から係数拡大の写像を作る -/

/-- ★係数制限は塔をなす。 -/
theorem isScalarTower_compHom (S T M : Type) [CommRing S] [CommRing T] [Algebra S T]
    [AddCommGroup M] [Module T M] :
    letI : Module S M := Module.compHom M (algebraMap S T)
    IsScalarTower S T M := by
  letI : Module S M := Module.compHom M (algebraMap S T)
  constructor
  intro a b x
  show (a • b) • x = (algebraMap S T a) • (b • x)
  rw [Algebra.smul_def, mul_smul]

/-- ★★★★★★**半線形写像 `p : M → N` から `T ⊗_S M →ₗ[T] N` を作る**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★これが `μ_L : 𝓞_K ⊗_{𝓞_F} Γ(L) → Γ(f^*L)` の本体である。 -/
noncomputable def bcLift (S T M N : Type) [CommRing S] [CommRing T] [Algebra S T]
    [AddCommGroup M] [Module S M] [AddCommGroup N] [Module T N]
    (p : M → N) (hadd : ∀ x y, p (x + y) = p x + p y)
    (hsmul : ∀ (a : S) (x : M), p (a • x) = (algebraMap S T a) • p x) :
    T ⊗[S] M →ₗ[T] N :=
  letI : Module S N := Module.compHom N (algebraMap S T)
  haveI : IsScalarTower S T N := isScalarTower_compHom S T N
  LinearMap.liftBaseChange T
    { toFun := p, map_add' := hadd, map_smul' := fun a x => hsmul a x }

@[simp] theorem bcLift_tmul (S T M N : Type) [CommRing S] [CommRing T] [Algebra S T]
    [AddCommGroup M] [Module S M] [AddCommGroup N] [Module T N]
    (p : M → N) (hadd : ∀ x y, p (x + y) = p x + p y)
    (hsmul : ∀ (a : S) (x : M), p (a • x) = (algebraMap S T a) • p x)
    (b : T) (x : M) :
    bcLift S T M N p hadd hsmul (b ⊗ₜ[S] x) = b • p x := by
  letI : Module S N := Module.compHom N (algebraMap S T)
  haveI : IsScalarTower S T N := isScalarTower_compHom S T N
  show LinearMap.liftBaseChange T
    ({ toFun := p, map_add' := hadd, map_smul' := fun a x => hsmul a x } : M →ₗ[S] N)
      (b ⊗ₜ[S] x) = _
  rw [LinearMap.liftBaseChange_tmul]
  rfl

/-! ### ★出典の紐付け(`.src`) -/

def pullSecTop_smul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(切断の引き戻しは係数環の側で半線形である)",
    sectionId := "genell-def-1-1-ii" }

def bcLift.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(半線形写像から係数拡大の写像 T ⊗_S M → N を作る)",
    sectionId := "genell-def-1-1-ii" }

def bcLift.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "pullSec_smul(切断の引き戻しは ρ-半線形、§9-741)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.pullSec_smul") 3,
    .citation "[mathlib]" "Scheme.ΓSpecIso_inv_naturality(Γ-Spec 同型の自然性)"
      (.inMathlib "AlgebraicGeometry.Scheme.ΓSpecIso_inv_naturality") 4,
    .citation "[mathlib]" "LinearMap.liftBaseChange(係数拡大の普遍性)"
      (.inMathlib "LinearMap.liftBaseChange") 4,
    .implicitStep
      ("★★これは底変換の**有限側**のための道具立てである。" ++
       "★★★残るのは μ_L が同型であること——道は §9-785 に書いたとおり " ++
       "μ_A ⊗ μ_{A⁻¹} ≅ μ_{Ō}(在庫の delta_pullSec と pullbackUnitPreIso)で回す") 4 ]

end ABC3.Found.Arakelov
