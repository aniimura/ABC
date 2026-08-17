import ABC3.Found.Arakelov.PicWhiskerW
import ABC3.Found.Arakelov.PicSheafTensor
import Mathlib.Algebra.Category.ModuleCat.Sheaf.Localization

/-!
# Arakelov (B1) 第 11 ブロック —— **層化が `f ▷ M` を同型に送る**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★`W` から `IsIso` へ

第 10 ブロックで `W (f ▷ M)` を取った。★★mathlib は

    (sheafification α).IsLocalization (J.W.inverseImage (toPresheaf R₀))

を持つ(`ModuleCat/Sheaf/Localization.lean`)ので、
`Localization.inverts` で **`IsIso` に変換できる**。

## ★★★★実装の大原則(この turn で 4 度踏んだ)

★★★**インスタンス束縛子は「型の書き方の違い」をまたげない。**

    X.PresheafOfModules      と  PresheafOfModules (X.presheaf ⋙ forget₂ _ _)

は `rfl` で等しいが、★**インスタンス探索は片方でしか成功しない**
(`Presheaf.IsLocallyInjective J (toSheafify J P.presheaf)` は前者でのみ見つかる)。

★★★**したがって局所単射性・局所全射性は `[..]` ではなく `(h : ..)` で受ける。**
そうすれば呼び出し側で形を合わせるだけで済む。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite

variable (X : Scheme.{u})

/-! ## ★★層化の単位は局所全単射 -/

/-- ★**層化の単位**(随伴の unit ではなく `toSheafify` として取る)。

★★★随伴の unit には局所単射・局所全射のインスタンスが付いていない(2026-08-17 実測)。 -/
noncomputable abbrev sheafifyUnit (P : X.PresheafOfModules) :
    P ⟶ ((PresheafOfModules.sheafification (R := X.ringCatSheaf)
      (𝟙 X.ringCatSheaf.obj)).obj P).val :=
  PresheafOfModules.toSheafify (𝟙 X.ringCatSheaf.obj)
    (CategoryTheory.toSheafify (Opens.grothendieckTopology X) P.presheaf)

/-- ★層化の単位は局所単射である。 -/
theorem isLocallyInjective_sheafifyUnit (P : X.PresheafOfModules) :
    PresheafOfModules.IsLocallyInjective (Opens.grothendieckTopology X)
      (PresheafOfModules.toSheafify (𝟙 X.ringCatSheaf.obj)
        (CategoryTheory.toSheafify (Opens.grothendieckTopology X) P.presheaf)) := by
  infer_instance

/-- ★層化の単位は局所全射である。 -/
theorem isLocallySurjective_sheafifyUnit (P : X.PresheafOfModules) :
    PresheafOfModules.IsLocallySurjective (Opens.grothendieckTopology X)
      (PresheafOfModules.toSheafify (𝟙 X.ringCatSheaf.obj)
        (CategoryTheory.toSheafify (Opens.grothendieckTopology X) P.presheaf)) := by
  infer_instance

/-! ## ★★可逆層の条件 -/

/-- ★**局所的に階数 1 自由**(可逆層の条件を、第 9 ブロックの仮定の形で述べたもの)。 -/
abbrev IsLocallyRankOne (M : X.PresheafOfModules) : Prop :=
  ∀ U : X.Opens, ∃ S : Sieve U, S ∈ (Opens.grothendieckTopology X) U ∧
    ∀ ⦃V : X.Opens⦄ (i : V ⟶ U), S i →
      Nonempty ((X.presheaf.obj (op V) : Type u) ≃ₗ[(X.presheaf.obj (op V) : Type u)]
        M.obj (op V))

/-! ## ★★★★`W` から `IsIso` へ -/

/-- ★★★★**層化は `f ▷ M` を同型に送る**。 -/
theorem isIso_sheafify_whiskerRight {P Q : X.PresheafOfModules} (f : P ⟶ Q)
    (hinj : PresheafOfModules.IsLocallyInjective (Opens.grothendieckTopology X) f)
    (hsur : PresheafOfModules.IsLocallySurjective (Opens.grothendieckTopology X) f)
    (M : X.PresheafOfModules) (hM : IsLocallyRankOne X M) :
    IsIso ((PresheafOfModules.sheafification (R := X.ringCatSheaf)
      (𝟙 X.ringCatSheaf.obj)).map (f ▷ M)) :=
  Localization.inverts _
    ((Opens.grothendieckTopology X).W.inverseImage
      (PresheafOfModules.toPresheaf X.ringCatSheaf.obj)) _
    (W_whiskerRight_of_basis (R := X.presheaf) (J := Opens.grothendieckTopology X)
      (f := f) (hinj := hinj) (hsur := hsur) (M := M) (htriv := hM))

/-- ★★★★**左からのテンソルでも同型**。 -/
theorem isIso_sheafify_whiskerLeft {P Q : X.PresheafOfModules} (f : P ⟶ Q)
    (hinj : PresheafOfModules.IsLocallyInjective (Opens.grothendieckTopology X) f)
    (hsur : PresheafOfModules.IsLocallySurjective (Opens.grothendieckTopology X) f)
    (M : X.PresheafOfModules) (hM : IsLocallyRankOne X M) :
    IsIso ((PresheafOfModules.sheafification (R := X.ringCatSheaf)
      (𝟙 X.ringCatSheaf.obj)).map (M ◁ f)) :=
  Localization.inverts _
    ((Opens.grothendieckTopology X).W.inverseImage
      (PresheafOfModules.toPresheaf X.ringCatSheaf.obj)) _
    (W_whiskerLeft_of_basis (R := X.presheaf) (J := Opens.grothendieckTopology X)
      (f := f) (hinj := hinj) (hsur := hsur) (M := M) (htriv := hM))

/-! ## ★★★★★可逆層とのテンソルは層化を変えない -/

/-- ★★★★**右から可逆層をテンソルしても層化は変わらない**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが結合律の 1 歩である。 -/
noncomputable def sheafifyTensorRight (P C : X.PresheafOfModules)
    (hC : IsLocallyRankOne X C) :
    (PresheafOfModules.sheafification (R := X.ringCatSheaf)
        (𝟙 X.ringCatSheaf.obj)).obj (P ⊗ C)
      ≅ (PresheafOfModules.sheafification (R := X.ringCatSheaf)
          (𝟙 X.ringCatSheaf.obj)).obj
          (((PresheafOfModules.sheafification (R := X.ringCatSheaf)
            (𝟙 X.ringCatSheaf.obj)).obj P).val ⊗ C) :=
  haveI := isIso_sheafify_whiskerRight X (sheafifyUnit X P)
    (isLocallyInjective_sheafifyUnit X P) (isLocallySurjective_sheafifyUnit X P) C hC
  asIso ((PresheafOfModules.sheafification (R := X.ringCatSheaf)
    (𝟙 X.ringCatSheaf.obj)).map ((sheafifyUnit X P) ▷ C))

/-- ★★★★**左から可逆層をテンソルしても層化は変わらない**。 -/
noncomputable def sheafifyTensorLeft (A P : X.PresheafOfModules)
    (hA : IsLocallyRankOne X A) :
    (PresheafOfModules.sheafification (R := X.ringCatSheaf)
        (𝟙 X.ringCatSheaf.obj)).obj (A ⊗ P)
      ≅ (PresheafOfModules.sheafification (R := X.ringCatSheaf)
          (𝟙 X.ringCatSheaf.obj)).obj
          (A ⊗ ((PresheafOfModules.sheafification (R := X.ringCatSheaf)
            (𝟙 X.ringCatSheaf.obj)).obj P).val) :=
  haveI := isIso_sheafify_whiskerLeft X (sheafifyUnit X P)
    (isLocallyInjective_sheafifyUnit X P) (isLocallySurjective_sheafifyUnit X P) A hA
  asIso ((PresheafOfModules.sheafification (R := X.ringCatSheaf)
    (𝟙 X.ringCatSheaf.obj)).map (A ◁ (sheafifyUnit X P)))

/-! ## ★出典の紐付け(`.src`) -/

def isIso_sheafify_whiskerRight.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——層化が f ▷ M を同型に送ること)",
    sectionId := "genell-def-1-1-i" }

def sheafifyTensorRight.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——可逆層とのテンソルは層化を変えないこと)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
