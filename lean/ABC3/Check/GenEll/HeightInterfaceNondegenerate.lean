import ABC3.Found.GenEll.HeightInterface

/-!
# `PulledBackClassData` の witness の非退化検査

`Found/GenEll/HeightInterface.lean` の docstring は
「本ファイルはその退化 witness ではない」と主張している。
★**散文のままでは検査されていない**ので、機械にかける。

**これは原典の主張ではない**(我々のモデルについての事実)ので `.src` を持たない。

## ★★★何を確かめるのか

G2 は「`PulledBackClassData.nonvacuous` という宣言があるか」しか見ていない。
★★`structure` の公理(`base_change_invariant` / `height_eq`)は
**`degOver` を恒等的に `0` にしても満たされる**——つまり
**G2 を通っただけでは高さが自明でないことは何も言えない**。

★★★本ファイルは 2 つを対にして測る:

1. **`structure` だけでは自明性を排除できない**(`trivialData`)——
   これが B5 の穴の正体である。
2. **我々の witness は自明ではない**(`heightData_degOver_eq`)——
   `degOver` は**原文の式そのもの** `deg_F(x_F^* M̄)` を返す。

★加えて、点の型が空でないこと(空なら `∀` がすべて空虚に真)も確かめる。
-/

namespace ABC3.Check.GenEll

open AlgebraicGeometry CategoryTheory NumberField
open ABC3.Interface.GenEll ABC3.Found.GenEll

/-! ## ★★★1. `structure` だけでは自明性を排除できない -/

/-- ★**すべてを `0` に送る `PulledBackClassData`**。

★★これが作れてしまうことが、`Interface` の docstring が言う
「退化 witness なら今すぐ作れる」の中身である。 -/
def trivialData : PulledBackClassData.{0, 0} where
  Point := PUnit
  Bundle := PUnit
  DefinedOver := fun _ _ _ _ => True
  degOver := fun _ _ _ _ _ => 0
  base_change_invariant := fun _ _ _ _ _ _ _ _ _ _ _ => rfl
  height := fun _ _ => 0
  height_eq := fun _ _ _ _ _ _ => rfl

/-- ★★★**`PulledBackClassData` の公理は高さの非自明性を強制しない**。

★これが `tools/check.mjs` 冒頭 B5 が名指しする穴である——
**G2 が通ったことは「中身がある」ことを意味しない**。 -/
theorem structure_alone_permits_zero :
    ∃ P : PulledBackClassData.{0, 0},
      ∀ (F : Type) [Field F] [NumberField F] (M : P.Bundle) (x : P.Point),
        P.degOver F M x = 0 :=
  ⟨trivialData, fun _ _ _ _ _ => rfl⟩

/-! ## ★★★2. 我々の witness は自明ではない -/

/-- ★★★**`degOver` は原文の式そのものを返す**。

原文 (GenEll p.4) の

    `ht_M̄(x) ≝ deg_F(x_F^* M̄)`

★★`F` 上のモデル `xF` を持つ点では、`degOver F` は
**`degNormalized (pullbackADiv F M̄ xF)`**、すなわち
`x_F^* M̄` の正規化次数**そのもの**である。

★★★**定数関数ではない**——これが `trivialData` との違いである。 -/
theorem heightData_degOver_eq {X : Scheme.{0}} (D₀ : X.IdealSheafData)
    (F : Type) [Field F] [NumberField F]
    (g : {g : GreenFn X // IsConjInvariant g})
    (xF : specRingOfIntegers F ⟶ X) (hJ : pullbackIdeal F D₀ xF ≠ 0) :
    (heightData D₀).degOver F g (UPoint.mk (algPointOff F xF hJ))
      = degNormalized (pullbackADiv F ⟨D₀, g.1⟩ xF) :=
  degOverField_eq ⟨D₀, g.1⟩ g.2 F _ (definedOverField_mk D₀ F xF hJ)

/-- ★★**`height` も同じ式である**(`degOver` と一致することは `height_eq`)。 -/
theorem heightData_height_eq {X : Scheme.{0}} (D₀ : X.IdealSheafData)
    (F : Type) [Field F] [NumberField F]
    (g : {g : GreenFn X // IsConjInvariant g})
    (xF : specRingOfIntegers F ⟶ X) (hJ : pullbackIdeal F D₀ xF ≠ 0) :
    (heightData D₀).height g (UPoint.mk (algPointOff F xF hJ))
      = degNormalized (pullbackADiv F ⟨D₀, g.1⟩ xF) :=
  rfl

/-! ## ★★3. 点の型が空でないこと -/

/-- ★★★**`Point` が空なら公理はすべて空虚に真になる**。

★供給した witness の `Point` は空でない——`Spec 𝓞_ℚ` の恒等射が点を与える。 -/
theorem heightData_point_nonempty :
    Nonempty (heightData (X := specRingOfIntegers ℚ)
      (⊤ : (specRingOfIntegers ℚ).IdealSheafData)).Point :=
  nonempty_uPoint_specRingOfIntegers

/-- ★★`DefinedOver` も空虚ではない——`ℚ` の上で定義される点が実際にある。

★★★これが無いと `base_change_invariant` と `height_eq` は
**前提が偽で自明に真**になってしまう。 -/
theorem heightData_definedOver_nonvacuous :
    ∃ x : (heightData (X := specRingOfIntegers ℚ)
      (⊤ : (specRingOfIntegers ℚ).IdealSheafData)).Point,
      (heightData (X := specRingOfIntegers ℚ)
        (⊤ : (specRingOfIntegers ℚ).IdealSheafData)).DefinedOver ℚ x := by
  refine ⟨UPoint.mk (algPointOff ℚ (𝟙 _) ?_), definedOverField_mk _ ℚ (𝟙 _) ?_⟩ <;>
    · rw [pullbackIdeal_top]
      exact fun h => absurd (h ▸ Submodule.mem_top (x := (1 : 𝓞 ℚ))) (by simp)

end ABC3.Check.GenEll
