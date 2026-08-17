import ABC3.Found.GenEll.ArchPoint

/-!
# [GenEll] Definition 1.2 / Proposition 1.4, (i) —— **高さを構成する**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5–6。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★高さを posit ではなく構成した

`Interface/GenEll/HeightTheory.lean` の `HeightTheoryData` は
**公理を 1 つも持たないデータ**であり、そのために
`Skeleton` の `prop_1_4` は**偽**になっていた
(反例は `Check/GenEll/HeightAxiomGap.lean`)。

★★本ファイルは `ht` を**実際に構成する**:

  `htArith D̄ x_F ≝ deg_F(x_F^* D̄)`

`x_F^* D̄` は `ArchPoint.lean` の `pullbackADiv` で、`deg_F` は `ArithDiv.lean` の
`degNormalized` である。★どちらも既に構成済みである。

## ★★★残る 1 本の等式

`Proposition 1.4, (i)`(加法性)に必要なものは、**ただ 1 本の等式**に落ちた:

  `x_F^*(D·E) = x_F^*D · x_F^*E`   —— `PullbackMul`

★★他はすべて証明済みである:
- アルキメデス側の加法性 —— `archADiv_add`(**無条件**)
- 次数の加法性 —— `deg_idealADiv_mul`(**無条件**)
- 因子の加法性 —— `idealADiv_mul`(**無条件**)

## ★★これは posit ではない

`PullbackMul` は `Interface` の `structure` ではなく、
**既に構成された対象についての等式**である。
★数学的には真である(`comap I f` は逆像イデアル層 `I·O_X` であり、
逆像イデアル層は積を保つ)。★★**mathlib に証明が無いだけである**——
`IdealSheaf/Functorial.lean` にあるのは `comap_sup`(左随伴なので上限は保つ)であり、
`comap_mul` は無い(2026-08-17 実測)。

★★★**したがって仮説として型に出し、非空虚性を示す**(`pullbackMul_id`)。
これは B5(条件を posit して未証明箇所を消す穴)ではない——
未証明箇所を隠すために足したのではなく、**残る 1 本を名指しした**のである。

## ★mathlib への貢献候補

`Scheme.IdealSheafData.comap_mul`。証明はアフィン局所で
`Ideal.map_mul` に落ちるが、mathlib の `comap` の定義
(ファイバー積の射影の核)からアフィン局所の記述を得る補題が無い。
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

/-- ★**空の算術因子の高さは 0**(非空虚性の witness)。 -/
@[simp] theorem htArith_one {X : Scheme.{0}} (xF : specRingOfIntegers F ⟶ X) :
    htArith F (ArithCartier.one X) xF = 0 := by
  rw [htArith, pullbackADiv_one, degNormalized_zero]

/-- ★★★**引き戻しの加法性** —— 残る 1 本を仮説に置いた形。

★★アルキメデス側は**無条件**、有限素点側は `PullbackMul` と非零性のみに依る。 -/
theorem pullbackADiv_tensor {X : Scheme.{0}} (D E : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X) (hmul : PullbackMul F xF)
    (hD : pullbackIdeal F D.divisor xF ≠ 0) (hE : pullbackIdeal F E.divisor xF ≠ 0) :
    pullbackADiv F (D.tensor E) xF = pullbackADiv F D xF + pullbackADiv F E xF := by
  refine Prod.ext ?_ ?_
  · show (pullbackADiv F (D.tensor E) xF).fin
      = (pullbackADiv F D xF + pullbackADiv F E xF).fin
    rw [ADiv.fin_add, pullbackADiv_fin, pullbackADiv_fin, pullbackADiv_fin]
    show (idealADiv F (pullbackIdeal F (D.divisor * E.divisor) xF)).fin = _
    rw [hmul D.divisor E.divisor, idealADiv_mul F _ _ hD hE, ADiv.fin_add]
  · show (pullbackADiv F (D.tensor E) xF).arc
      = (pullbackADiv F D xF + pullbackADiv F E xF).arc
    rw [ADiv.arc_add]
    exact pullbackADiv_arc_tensor F D E xF

/-- ★★★**`Proposition 1.4, (i)`** —— 構成した高さの加法性。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

★★★**posit された `HeightTheoryData` の上ではこれは偽だった**
(`Check/GenEll/HeightAxiomGap.lean` に反例)。
構成に置き換えると**真になる**——それが本ファイルの主眼である。 -/
theorem htArith_tensor {X : Scheme.{0}} (D E : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X) (hmul : PullbackMul F xF)
    (hD : pullbackIdeal F D.divisor xF ≠ 0) (hE : pullbackIdeal F E.divisor xF ≠ 0) :
    htArith F (D.tensor E) xF = htArith F D xF + htArith F E xF := by
  rw [htArith, htArith, htArith, pullbackADiv_tensor F D E xF hmul hD hE,
    degNormalized_add]

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

def htArith_tensor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (i)(構成した高さの加法性——引き戻しの積保存を仮説に置く)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
