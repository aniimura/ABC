/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.BaseChangeUnit
import ABC3.Found.Arakelov.SpanPullSecInv
import ABC3.Found.Arakelov.GammaModInvertible
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★**底変換の同型** `𝓞_K ⊗_{𝓞_F} Γ(L) ≅ Γ(f^*L)`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

## ★★★★★★★★★★これは何か

    `μ_L : T ⊗_S Γ(L) ⟶ Γ(f^*L)`   は**同型**である（`L` は算術直線束）。

★これで `deg_K(L|_{Spec 𝓞_K}) = deg_F(L)` の**有限側**が計算できるようになる。

## ★★★★★機構（3 段、すべて本セッションで積んだ）

    span_pullSecT_of_inv          : `Γ(f^*L)` は引き戻した大域切断で生成される（`§9-790`）
    invertible_gammaModPre        : `Γ(L)` は係数環の上でも可逆（`§9-788`）
    Module.Invertible.bijective_of_surjective : 可逆加群の間の全射は全単射（mathlib）

★★`§9-790` は `Γ(X,⊤)`-span で述べられているので、
`Submodule.span_induction` で `μ` の像に落とす。
★★★係数の食い違い（`T`-作用 vs `Γ(X,⊤)`-作用）は `gammaModPre_smul`（`rfl`）で吸収する
——`gammaModPre` を `abbrev` にしたのでインスタンス探索が通る。

## ★★★これで何ができるか

    `Γ(f^*L) / Γ_X·(pullSec s) ≅ T ⊗_S (Γ(L) / Γ_Y·s)`

が言えるので、`#(T ⊗_S Q) = #(Q)^{[K:F]}` と合わせれば底変換の有限側が出る。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace
open ABC3.Found.GenEll
open scoped TensorProduct

/-- ★`gammaModPre` の係数作用は `Γ-Spec` 同型で書ける（`rfl`）。 -/
theorem gammaModPre_smul (R : CommRingCat.{0}) (M : (Spec R).PresheafOfModules)
    (t : (R : Type)) (x : (gammaModPre R M : Type)) :
    t • x = ((Scheme.ΓSpecIso R).inv.hom t) • x := rfl

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★**`μ_L` は全射である**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★`§9-790` の生成性（`Γ(X,⊤)`-span）を `Submodule.span_induction` で
`μ` の像に落とすだけである。 -/
theorem muBC_surjective (S T : Type) [CommRing S] [CommRing T] [Algebra S T]
    (L : AInv (Spec (CommRingCat.of S))) :
    Function.Surjective (muBC S T L.carrier.sheaf) := by
  intro y
  have hy : y ∈ Submodule.span (Γ(Spec (CommRingCat.of T),
        (⊤ : (Spec (CommRingCat.of T)).Opens)) : Type)
      (Set.range (pullSecT (Spec.map (CommRingCat.ofHom (algebraMap S T))) L.carrier.sheaf)) := by
    rw [span_pullSecT_of_inv]; trivial
  have hgoal : y ∈ LinearMap.range (muBC S T L.carrier.sheaf) := by
    induction hy using Submodule.span_induction with
    | mem x hx =>
        obtain ⟨s, rfl⟩ := hx
        exact ⟨(1 : T) ⊗ₜ[S] s, by
          rw [muBC, bcLift_tmul, one_smul]
          rfl⟩
    | zero => exact Submodule.zero_mem _
    | add u v hu hv hu' hv' => exact Submodule.add_mem _ hu' hv'
    | smul c x hx hx' =>
        obtain ⟨z, hz⟩ := hx'
        refine ⟨((Scheme.ΓSpecIso (CommRingCat.of T)).hom.hom c) • z, ?_⟩
        rw [map_smul, hz, gammaModPre_smul]
        have hc : (Scheme.ΓSpecIso (CommRingCat.of T)).inv.hom
            ((Scheme.ΓSpecIso (CommRingCat.of T)).hom.hom c) = c :=
          congrArg (fun (m : _ ⟶ _) => CommRingCat.Hom.hom m c)
            (Scheme.ΓSpecIso (CommRingCat.of T)).hom_inv_id
        rw [hc]
        rfl
  exact hgoal

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★★**`μ_L` は全単射である**——底変換の同型。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★両側とも可逆加群（`§9-788`）なので、全射から全単射に上がる。 -/
theorem muBC_bijective (S T : Type) [CommRing S] [CommRing T] [Algebra S T]
    (L : AInv (Spec (CommRingCat.of S))) :
    Function.Bijective (muBC S T L.carrier.sheaf) := by
  haveI h1 : Module.Invertible S ((gammaModPre (CommRingCat.of S) L.carrier.sheaf) : Type) :=
    invertible_gammaModPre (CommRingCat.of S) L
  haveI h2 : Module.Invertible T ((gammaModPre (CommRingCat.of T)
      ((pullbackPre (Spec.map (CommRingCat.ofHom (algebraMap S T)))).obj
        L.carrier.sheaf)) : Type) :=
    invertible_gammaModPre (CommRingCat.of T)
      (AInv.pullback (Spec.map (CommRingCat.ofHom (algebraMap S T))) L)
  exact Module.Invertible.bijective_of_surjective (muBC_surjective S T L)

/-- ★★★★★★★★★★**底変換の同型** `T ⊗_S Γ(L) ≃ₗ[T] Γ(f^*L)`。 -/
noncomputable def muBCEquiv (S T : Type) [CommRing S] [CommRing T] [Algebra S T]
    (L : AInv (Spec (CommRingCat.of S))) :
    T ⊗[S] ((gammaModPre (CommRingCat.of S) L.carrier.sheaf) : Type)
      ≃ₗ[T] ((gammaModPre (CommRingCat.of T)
        ((pullbackPre (Spec.map (CommRingCat.ofHom (algebraMap S T)))).obj
          L.carrier.sheaf)) : Type) :=
  LinearEquiv.ofBijective (muBC S T L.carrier.sheaf) (muBC_bijective S T L)

/-- ★切断の引き戻し（係数環の側で見た型で宣言し直したもの）。 -/
noncomputable def pullSecMod (S T : Type) [CommRing S] [CommRing T] [Algebra S T]
    (L : (Spec (CommRingCat.of S)).PresheafOfModules)
    (x : ((gammaModPre (CommRingCat.of S) L) : Type)) :
    ((gammaModPre (CommRingCat.of T)
      ((pullbackPre (Spec.map (CommRingCat.ofHom (algebraMap S T)))).obj L)) : Type) :=
  pullSecTop (CommRingCat.ofHom (algebraMap S T)) L x

theorem muBCEquiv_tmul (S T : Type) [CommRing S] [CommRing T] [Algebra S T]
    (L : AInv (Spec (CommRingCat.of S))) (b : T)
    (x : ((gammaModPre (CommRingCat.of S) L.carrier.sheaf) : Type)) :
    muBCEquiv S T L (b ⊗ₜ[S] x) = b • (pullSecMod S T L.carrier.sheaf x) :=
  bcLift_tmul S T _ _ _ _ _ b x

/-! ### ★出典の紐付け(`.src`) -/

def muBC_bijective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(底変換の同型 𝓞_K ⊗_{𝓞_F} Γ(L) ≅ Γ(f^*L))",
    sectionId := "genell-def-1-1-ii" }

def muBC_bijective.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "span_pullSecT_of_inv(Γ(f^*L) は引き戻した大域切断で生成される、§9-790)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.span_pullSecT_of_inv") 3,
    .citation "[ABC3]" "invertible_gammaModPre(Γ(L) は係数環の上でも可逆、§9-788)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.invertible_gammaModPre") 4,
    .citation "[mathlib]" "Module.Invertible.bijective_of_surjective"
      (.inMathlib "Module.Invertible.bijective_of_surjective") 4,
    .implicitStep
      ("★★これで deg_K(L|_{Spec 𝓞_K}) = deg_F(L) の**有限側**が計算できる: " ++
       "Γ(f^*L) / Γ_X·(pullSec s) ≅ T ⊗_S (Γ(L) / Γ_Y·s) となるので、" ++
       "#(T ⊗_S Q) = #(Q)^{[K:F]} と合わせればよい") 4 ]

end ABC3.Found.Arakelov
