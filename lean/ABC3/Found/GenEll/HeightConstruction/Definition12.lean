/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ArchPoint

/-!
# HeightConstruction —— `[GenEll] Definition 1.2` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField
variable (F : Type) [Field F] [NumberField F]

/-! ## ★★残る 1 本を型に出す -/

/-- ★★★**引き戻しが積を保つ**という条件。

★これは `Proposition 1.4, (i)` に残った**ただ 1 本**である。
数学的には真だが、mathlib に `IdealSheafData.comap_mul` が無い。 -/
def PullbackMul {X : Scheme.{0}} (xF : specRingOfIntegers F ⟶ X) : Prop :=
  ∀ D E : X.IdealSheafData,
    pullbackIdeal F (D * E) xF = pullbackIdeal F D xF * pullbackIdeal F E xF

/-- ★★**非空虚性の witness** —— 恒等射では `PullbackMul` が**証明できる**。

★機構は 3 つ:
`comap_id`(恒等射の引き戻しは何もしない)、
`equivOfIsAffine` が**環同型**であること(`≃+*o` なので `map_mul` が効く)、
`Ideal.comap_symm`(同型に沿った `comap` は `map` に書き換わり、`Ideal.map_mul` が効く)。

★★**証明は `NumberField F` を使っていない**(linter が検出)——
つまり**任意のアフィンスキームで成り立つ**。恒等射に限った現象ではなく、
`comap` が問題になるのは**射が非自明なとき**だけである。 -/
theorem pullbackMul_id : PullbackMul F (𝟙 (specRingOfIntegers F)) := by
  intro D E
  have hcomap : ∀ I : (specRingOfIntegers F).IdealSheafData,
      pullbackIdeal F I (𝟙 (specRingOfIntegers F))
        = (Scheme.IdealSheafData.equivOfIsAffine I).comap
            (Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).inv.hom := by
    intro I
    rw [pullbackIdeal, Scheme.IdealSheafData.comap_id]
  -- `comap` を同型の逆向きの `map` に書き換える
  set e : Γ(specRingOfIntegers F, ⊤) ≃+* (𝓞 F) :=
    (Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).commRingCatIsoToRingEquiv with he
  -- ★同型に沿った `comap` は逆向きの `map` である(定義的に等しいので `rfl` で渡る)
  have key : ∀ J : Ideal Γ(specRingOfIntegers F, ⊤),
      Ideal.comap (Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).inv.hom J
        = J.map (e : Γ(specRingOfIntegers F, ⊤) →+* (𝓞 F)) :=
    fun J => Ideal.comap_symm e
  rw [hcomap, hcomap, hcomap, key, key, key, map_mul, Ideal.map_mul]

/-! ## ★★★高さの構成 -/

/-- ★★★**高さ** `ht_{D̄}(x) ≝ deg_F(x_F^* D̄)`。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★★**posit ではなく構成である。**
`pullbackADiv`(`ArchPoint.lean`)と `degNormalized`(`ArithDiv.lean`)の合成。 -/
noncomputable def htArith {X : Scheme.{0}} (D : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X) : ℝ :=
  degNormalized (pullbackADiv F D xF)

/-- ★★**アルキメデス側だけを見れば加法性は無条件に成り立つ**。

★仮説がどこに効いているかを型で分けたもの。 -/
theorem htArith_arc_unconditional {X : Scheme.{0}} (D E : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X) :
    (pullbackADiv F (D.tensor E) xF).arc
      = (pullbackADiv F D xF).arc + (pullbackADiv F E xF).arc :=
  pullbackADiv_arc_tensor F D E xF

/-! ## ★出典の紐付け(`.src`)

★条つきである。`Definition 1.2` 全体には `Definition 1.1`(算術直線束そのもの)が要り、
`Proposition 1.4` 全体には (i)〜(iv) の 4 条すべてが要る。 -/

def htArith.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2, (i)(高さの構成——因子とグリーン関数による)",
    sectionId := "genell-def-1-2-i" }

end ABC3.Found.GenEll
