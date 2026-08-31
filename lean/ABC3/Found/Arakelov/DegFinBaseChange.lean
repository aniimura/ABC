/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.DegFinBCTools
import ABC3.Found.Arakelov.DegArithPre
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★底変換の**有限側が閉じた** `degFin_K(f^*L̄) = [K:F]·degFin_F(L̄)`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> degK(L|Spec(OK)) = degF (L)

## ★★★★★★★★★★これは何か

    `degFin_K(f^*L̄, f^*s) = [K:F] · degFin_F(L̄, s)`

★正規化 `1/[K:ℚ] = 1/([F:ℚ]·[K:F])` がこの `[K:F]` を打ち消す
——これが原文の `deg_K(L|_{Spec 𝓞_K}) = deg_F(L)` の有限部分である。

## ★★★★★機構（5 段、すべて本セッションで積んだ）

    card_quotient_span_gammaModPre : `Γ(X,⊤)`-span と `R`-span は同じ（`§9-797`）
    muBCEquiv                      : `Γ(f^*L) ≅ 𝓞_K ⊗_{𝓞_F} Γ(L)`（`§9-791`）
    card_quotient_tensor_span      : `#((T⊗A)/T·(1⊗a)) = #(T⊗(A/R·a))`（`§9-797`）
    card_tensor_ringOfIntegers     : `#(𝓞_K ⊗ Q) = #(Q)^{[K:F]}`（`§9-794`）
    Real.log_pow                   : `log (x^n) = n·log x`

★★`muBCEquiv` が `1 ⊗ s` を `pullSec s` に送る（`muBCEquiv_tmul`）ので、
商の対応がそのまま取れる。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace NumberField
open ABC3.Found.GenEll
open scoped TensorProduct

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★**商の位数は `[K:F]` 乗になる**。

原文 (GenEll p.4):
> degK(L|Spec(OK)) = degF (L) -/
theorem card_quotient_baseChange (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K]
    (L : AInv (Spec (CommRingCat.of (𝓞 F))))
    (s : (gammaModPre (CommRingCat.of (𝓞 F)) L.carrier.sheaf : Type)) (hs : s ≠ 0) :
    Nat.card ((gammaModPre (CommRingCat.of (𝓞 K))
        ((pullbackPre (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))))).obj
          L.carrier.sheaf) : Type)
        ⧸ ((Γ(Spec (CommRingCat.of (𝓞 K)), (⊤ : (Spec (CommRingCat.of (𝓞 K))).Opens)) : Type)
          ∙ (pullSecMod (𝓞 F) (𝓞 K) L.carrier.sheaf s)))
      = (Nat.card ((gammaModPre (CommRingCat.of (𝓞 F)) L.carrier.sheaf : Type)
          ⧸ ((Γ(Spec (CommRingCat.of (𝓞 F)), (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type)
            ∙ s))) ^ (Module.finrank F K) := by
  haveI hinvF := invertible_gammaModPre (CommRingCat.of (𝓞 F)) L
  haveI hQfin : Finite ((gammaModPre (CommRingCat.of (𝓞 F)) L.carrier.sheaf : Type)
      ⧸ (((CommRingCat.of (𝓞 F)) : Type) ∙ s)) :=
    finite_quotient_span_invertible ((CommRingCat.of (𝓞 F)) : Type)
      (gammaModPre (CommRingCat.of (𝓞 F)) L.carrier.sheaf : Type)
      (fun r hr => finite_quotient_span F r hr) s hs
  have hmap : Submodule.map ((muBCEquiv (𝓞 F) (𝓞 K) L) : _ →ₗ[(𝓞 K)] _)
      (Submodule.span (𝓞 K) {(1 : (𝓞 K)) ⊗ₜ[(𝓞 F)] s})
      = Submodule.span (𝓞 K) {pullSecMod (𝓞 F) (𝓞 K) L.carrier.sheaf s} := by
    rw [Submodule.map_span]
    simp only [Set.image_singleton, LinearEquiv.coe_coe]
    congr 1
    rw [muBCEquiv_tmul, one_smul]
  calc Nat.card (_ ⧸ ((Γ(Spec (CommRingCat.of (𝓞 K)),
        (⊤ : (Spec (CommRingCat.of (𝓞 K))).Opens)) : Type)
          ∙ (pullSecMod (𝓞 F) (𝓞 K) L.carrier.sheaf s)))
      = Nat.card (_ ⧸ (((CommRingCat.of (𝓞 K)) : Type)
          ∙ (pullSecMod (𝓞 F) (𝓞 K) L.carrier.sheaf s))) :=
        card_quotient_span_gammaModPre (CommRingCat.of (𝓞 K)) _ _
    _ = Nat.card (((𝓞 K) ⊗[(𝓞 F)]
          (gammaModPre (CommRingCat.of (𝓞 F)) L.carrier.sheaf : Type))
          ⧸ (Submodule.span (𝓞 K) {(1 : (𝓞 K)) ⊗ₜ[(𝓞 F)] s})) :=
        (Nat.card_congr (Submodule.Quotient.equiv _ _ (muBCEquiv (𝓞 F) (𝓞 K) L) hmap).toEquiv).symm
    _ = Nat.card ((𝓞 K) ⊗[(𝓞 F)]
          ((gammaModPre (CommRingCat.of (𝓞 F)) L.carrier.sheaf : Type)
            ⧸ (Submodule.span (𝓞 F) {s}))) :=
        card_quotient_tensor_span (𝓞 F) (𝓞 K) _ s
    _ = (Nat.card ((gammaModPre (CommRingCat.of (𝓞 F)) L.carrier.sheaf : Type)
          ⧸ (Submodule.span (𝓞 F) {s}))) ^ (Module.finrank F K) :=
        card_tensor_ringOfIntegers F K _
    _ = _ := by rw [card_quotient_span_gammaModPre (CommRingCat.of (𝓞 F)) L.carrier.sheaf s]

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★**底変換の有限側**——`degFin_K(f^*L̄) = [K:F]·degFin_F(L̄)`。

原文 (GenEll p.4):
> degK(L|Spec(OK)) = degF (L)

★正規化 `1/[K:ℚ] = 1/([F:ℚ]·[K:F])` がこの `[K:F]` を打ち消す。 -/
theorem degFinPre_baseChange (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K]
    (L : AInv (Spec (CommRingCat.of (𝓞 F))))
    (s : (gammaModPre (CommRingCat.of (𝓞 F)) L.carrier.sheaf : Type)) (hs : s ≠ 0) :
    degFinPre (AInv.pullback (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))) L)
        (pullSecMod (𝓞 F) (𝓞 K) L.carrier.sheaf s)
      = (Module.finrank F K : ℝ) * degFinPre L s := by
  show Real.log (Nat.card (_ ⧸ ((Γ(Spec (CommRingCat.of (𝓞 K)),
      (⊤ : (Spec (CommRingCat.of (𝓞 K))).Opens)) : Type)
        ∙ (pullSecMod (𝓞 F) (𝓞 K) L.carrier.sheaf s)))) = _
  rw [card_quotient_baseChange F K L s hs]
  push_cast
  rw [Real.log_pow]
  rfl

/-! ### ★出典の紐付け(`.src`) -/

def degFinPre_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(底変換の有限側——degFin_K = [K:F]·degFin_F)",
    sectionId := "genell-def-1-1-ii" }

def degFinPre_baseChange.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "muBCEquiv(Γ(f^*L) ≅ 𝓞_K ⊗_{𝓞_F} Γ(L)、§9-791)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.muBCEquiv") 4,
    .citation "[ABC3]" "card_tensor_ringOfIntegers(#(𝓞_K ⊗ Q) = #(Q)^{[K:F]}、§9-794)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.card_tensor_ringOfIntegers") 4,
    .citation "[ABC3]" "card_quotient_tensor_span / card_quotient_span_gammaModPre(§9-797)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.card_quotient_tensor_span") 4,
    .implicitStep
      ("★★これと §9-796 の archDeg_baseChange を足すと degArithPre の底変換不変性が出る。" ++
       "★★★そのとき pullSec s ≠ 0 が要る——𝓞_K が 𝓞_F 上忠実平坦であることから " ++
       "1 ⊗ s ≠ 0 が出るはずである(未着手)") 4 ]

end ABC3.Found.Arakelov
