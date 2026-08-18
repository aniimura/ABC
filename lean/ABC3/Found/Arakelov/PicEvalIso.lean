import ABC3.Found.Arakelov.PicEvalBij
import ABC3.Found.Arakelov.PicLocalBij
import ABC3.Found.Arakelov.PicSheafGroup

/-!
# Arakelov (B2) 第 182 ブロック —— ★★★★★★★**局所自明な層加群は可逆層である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★B2 の本丸が落ちた

    IsLocallyTrivial X F.val  ⟹  InvSheaf X

これは B1 の `invSheafOfModule`(第 133、`Spec R` 限定)の**一般スキーム版**であり、
`ofDivisor` を作るための土台である。

## ★★機構は 3 段

| 段 | 内容 | 在庫 |
|---|---|---|
| 1 | 評価射は自明な所で全単射 | ★第 181 |
| 2 | ゆえに**局所全単射** | 第 103(`isLocallyBijective_of_cover`) |
| 3 | 局所全単射は層化で同型 | ★mathlib `inverseImage_W_toPresheaf_eq_inverseImage_isomorphisms` |

★★★第 3 段は **2026-08-19 に実測**して見つけた:
`Mathlib/Algebra/Category/ModuleCat/Sheaf/Localization.lean` に
「層化は `J.W`(局所全単射)を同型に送る」が**morphism property の等式**として在る。
★これが無ければ層化の普遍性から手で組む羽目になっていた(見積もり +6 ブロック)。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `locallyBijective_evHom` | ★★評価射は局所全単射 |
| `isIso_sheafify_evHom` | ★★★★評価射の層化は同型 |
| `evalModulesIso` | ★★★★★`F ⊗ F^∨ ≅ 𝒪_X` |
| `InvSheaf.ofLocallyTrivial` | ★★★★★★★**局所自明な層加群は可逆層** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} (F : X.PresheafOfModules) (hF : IsLocallyTrivial X F)

include hF in
/-- ★★**評価射は局所全単射である**。 -/
theorem locallyBijective_evHom :
    Presheaf.IsLocallyInjective (Opens.grothendieckTopology X)
        ((PresheafOfModules.toPresheaf _).map (evHom F)) ∧
      Presheaf.IsLocallySurjective (Opens.grothendieckTopology X)
        ((PresheafOfModules.toPresheaf _).map (evHom F)) := by
  have hb : ∀ U : X.Opens, ∃ S : Sieve U, S ∈ (Opens.grothendieckTopology X) U ∧
      ∀ ⦃W : X.Opens⦄ (i : W ⟶ U), S.arrows i →
        Function.Bijective (((PresheafOfModules.toPresheaf _).map (evHom F)).app (op W)) := by
    intro U
    obtain ⟨S, hS, hiso⟩ := hF U
    refine ⟨S, hS, fun {W} i hi => ?_⟩
    obtain ⟨e⟩ := hiso i hi
    exact bijective_evHom_app F W e
  constructor
  · refine isLocallyInjective_of_coverSieve _ _ (fun U x y hxy => ?_)
    obtain ⟨S, hS, hbb⟩ := hb U
    refine ⟨S, hS, fun {W} i hi => ?_⟩
    refine (hbb i hi).1 ?_
    rw [NatTrans.naturality_apply ((PresheafOfModules.toPresheaf _).map (evHom F)) i.op x,
      NatTrans.naturality_apply ((PresheafOfModules.toPresheaf _).map (evHom F)) i.op y, hxy]
  · refine isLocallySurjective_of_cover _ _ (fun U s => ?_)
    obtain ⟨S, hS, hbb⟩ := hb U
    exact ⟨S, hS, fun {W} i hi => (hbb i hi).2 _⟩

include hF in
/-- ★★★★**評価射の層化は同型である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★mathlib の「層化は局所全単射を同型に送る」(`Sheaf/Localization.lean`)に乗せる。 -/
theorem isIso_sheafify_evHom : IsIso ((sheafifyFunctor X).map (evHom F)) := by
  obtain ⟨hi, hs⟩ := locallyBijective_evHom F hF
  haveI hi' : Presheaf.IsLocallyInjective (Opens.grothendieckTopology X)
      ((PresheafOfModules.toPresheaf X.ringCatSheaf.obj).map (evHom F)) := hi
  haveI hs' : Presheaf.IsLocallySurjective (Opens.grothendieckTopology X)
      ((PresheafOfModules.toPresheaf X.ringCatSheaf.obj).map (evHom F)) := hs
  have hW : (MorphismProperty.inverseImage (Opens.grothendieckTopology X).W
      (PresheafOfModules.toPresheaf X.ringCatSheaf.obj)) (evHom F) :=
    GrothendieckTopology.W_of_isLocallyBijective _ _
  have heq := PresheafOfModules.inverseImage_W_toPresheaf_eq_inverseImage_isomorphisms
    (R := X.ringCatSheaf) (𝟙 X.ringCatSheaf.obj)
  have h2 : (MorphismProperty.isomorphisms (SheafOfModules X.ringCatSheaf)).inverseImage
      (PresheafOfModules.sheafification (𝟙 X.ringCatSheaf.obj)) (evHom F) := heq ▸ hW
  exact h2

/-- ★★★★★**`F ⊗ F^∨ ≅ 𝒪_X`**。 -/
noncomputable def evalModulesIso (F : X.Modules) (hF : IsLocallyTrivial X F.val) :
    tensorModules F (dualModules F) ≅ unitModules X :=
  @asIso _ _ _ _ ((sheafifyFunctor X).map (evHom F.val)) (isIso_sheafify_evHom F.val hF)
    ≪≫ sheafifyValIso (unitModules X)

/-- ★★★★★★★**局所自明な層加群は可逆層である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★★第 133 の `invSheafOfModule`(`Spec R` 限定)の**一般スキーム版**である。
これが `ofDivisor`(因子から直線束)の土台になる。 -/
noncomputable def InvSheaf.ofLocallyTrivial (F : X.Modules) (hF : IsLocallyTrivial X F.val) :
    InvSheaf X where
  carrier := F
  inv := dualModules F
  isInv := ⟨evalModulesIso F hF⟩
  trivial := hF
  invTrivial := isLocallyTrivial_dualPresheaf F.val hF

/-! ## ★出典の紐付け(`.src`) -/

def InvSheaf.ofLocallyTrivial.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——局所自明な層加群は可逆層)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
