import ABC3.Found.Arakelov.PicLocalBasis

/-!
# Arakelov (B1) 第 10 ブロック —— **`W` は可逆層とのテンソルで保たれる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★これが結合律の中核である

層加群のテンソル積の結合律は

    層化 (層化 (A ⊗ B).val ⊗ C) ≅ 層化 ((A ⊗ B) ⊗ C)

に帰着し、これは層化の単位 `η : A ⊗ B ⟶ 層化(A ⊗ B).val` について

    W (η ▷ C)          (W = 層化が反転させる射のクラス)

を言うことである。

★★★`WEqualsLocallyBijective` により `W f ⟺ 局所単射 ∧ 局所全射`。
2 条はそれぞれ第 6・第 9 ブロックで取ってある。**本ブロックはその合流である。**
-/

universe u

namespace ABC3.Found.Arakelov

open CategoryTheory MonoidalCategory Opposite

/-- ★★★★★**`W` は可逆層とのテンソルで保たれる**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで結合律の中核が取れた。

★機構は第 6 ブロック(局所全射)+ 第 9 ブロック(局所単射)+
`GrothendieckTopology.W_of_isLocallyBijective`。 -/
theorem W_whiskerRight_of_basis {C : Type u} [Category.{u} C]
    (J : GrothendieckTopology C) [J.WEqualsLocallyBijective AddCommGrpCat.{u}]
    {R : Cᵒᵖ ⥤ CommRingCat.{u}}
    {P Q : PresheafOfModules.{u} (R ⋙ forget₂ CommRingCat.{u} RingCat.{u})}
    (f : P ⟶ Q)
    (hinj : PresheafOfModules.IsLocallyInjective J f)
    (hsur : PresheafOfModules.IsLocallySurjective J f)
    (M : PresheafOfModules.{u} (R ⋙ forget₂ CommRingCat.{u} RingCat.{u}))
    (htriv : ∀ U : C, ∃ S : Sieve U, S ∈ J U ∧ ∀ ⦃V : C⦄ (i : V ⟶ U), S i →
      Nonempty ((R.obj (op V) : Type u) ≃ₗ[(R.obj (op V) : Type u)] M.obj (op V))) :
    (J.W.inverseImage
        (PresheafOfModules.toPresheaf (R ⋙ forget₂ CommRingCat.{u} RingCat.{u}))) (f ▷ M) := by
  haveI := hinj
  haveI := hsur
  rw [MorphismProperty.inverseImage_iff]
  haveI := isLocallyInjective_whiskerRight_of_basis J f M htriv
  haveI := isLocallySurjective_whiskerRight J f M
  exact J.W_of_isLocallyBijective _

/-- ★★★**左からのテンソルでも `W` は保たれる**(対称性)。

★機構は mathlib の `whiskerRight` と同じ——`MorphismProperty.arrow_mk_iso_iff` に
braiding で作った射の同型を当てる。 -/
theorem W_whiskerLeft_of_basis {C : Type u} [Category.{u} C]
    (J : GrothendieckTopology C) [J.WEqualsLocallyBijective AddCommGrpCat.{u}]
    {R : Cᵒᵖ ⥤ CommRingCat.{u}}
    {P Q : PresheafOfModules.{u} (R ⋙ forget₂ CommRingCat.{u} RingCat.{u})}
    (f : P ⟶ Q)
    (hinj : PresheafOfModules.IsLocallyInjective J f)
    (hsur : PresheafOfModules.IsLocallySurjective J f)
    (M : PresheafOfModules.{u} (R ⋙ forget₂ CommRingCat.{u} RingCat.{u}))
    (htriv : ∀ U : C, ∃ S : Sieve U, S ∈ J U ∧ ∀ ⦃V : C⦄ (i : V ⟶ U), S i →
      Nonempty ((R.obj (op V) : Type u) ≃ₗ[(R.obj (op V) : Type u)] M.obj (op V))) :
    (J.W.inverseImage
        (PresheafOfModules.toPresheaf (R ⋙ forget₂ CommRingCat.{u} RingCat.{u}))) (M ◁ f) :=
  ((J.W.inverseImage (PresheafOfModules.toPresheaf _)).arrow_mk_iso_iff
      (Arrow.isoMk (β_ M P) (β_ M Q)
        (BraidedCategory.braiding_naturality_right M f).symm)).2
    (W_whiskerRight_of_basis J f hinj hsur M htriv)

/-! ## ★出典の紐付け(`.src`) -/

def W_whiskerRight_of_basis.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——W が可逆層とのテンソルで保たれること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
