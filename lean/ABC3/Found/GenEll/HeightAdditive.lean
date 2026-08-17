import ABC3.Found.GenEll.ComapMul
import ABC3.Found.GenEll.HeightConstruction

/-!
# [GenEll] Proposition 1.4, (i) —— **高さの加法性(無条件)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

## ★★★仮説が消えた

`HeightConstruction.lean` は `ht` の加法性を **`PullbackMul` 1 本**の仮説に落とした:

    `x_F^*(D·E) = x_F^*D · x_F^*E`

★その正体は `Scheme.IdealSheafData.comap_mul` であり、mathlib に無かった。
★★★**`ComapMul.lean` で証明したので、仮説は消える。**

## ★★残るのは原文自身の仮定だけ

`x_F^* D ≠ 0` は原文の「`x_F` は因子 `D` を通らない」に当たる
(`StalkPullback.lean` の `isEffectiveCartierAt_stalkPullback_of_ne_bot` を参照——
茎が整域なので `≠ ⊥` だけでよい)。

★★★**したがって本ファイルの `htArith_tensor_unconditional` は
`Proposition 1.4, (i)` そのものである。**

## ★★これで `Proposition 1.4` の 3 条が無条件になった

| 条 | 状態 |
|---|---|
| (i) 加法性 | ★★★**本ファイルで無条件に** |
| (ii) 非負性 | ★★無条件(`HeightNonneg.lean`) |
| (iii) 類にしか依らない | ★★無条件(`HeightClass.lean`) |
| (iv) Northcott | 未着手 |
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

variable (F : Type) [Field F] [NumberField F]

/-! ## ★★★`PullbackMul` は常に成り立つ -/

/-- ★★★**引き戻しは常に積を保つ**——仮説ではなく定理である。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

★機構は 3 段:
`comap_mul`(`ComapMul.lean`、本日取得)、
`equivOfIsAffine` が**環同型**であること、
`Ideal.comap_symm`(同型に沿った `comap` は `map` に書き換わる)。

★★`HeightConstruction.lean` の `pullbackMul_id`(恒等射の場合)を
**任意の射に一般化したもの**である。 -/
theorem pullbackMul_all {X : Scheme.{0}} (xF : specRingOfIntegers F ⟶ X) :
    PullbackMul F xF := by
  intro D E
  have hcomap : ∀ I : X.IdealSheafData,
      pullbackIdeal F I xF
        = (Scheme.IdealSheafData.equivOfIsAffine (I.comap xF)).comap
            (Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).inv.hom := fun _ => rfl
  set e : Γ(specRingOfIntegers F, ⊤) ≃+* (𝓞 F) :=
    (Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).commRingCatIsoToRingEquiv with he
  have key : ∀ J : Ideal Γ(specRingOfIntegers F, ⊤),
      Ideal.comap (Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).inv.hom J
        = J.map (e : Γ(specRingOfIntegers F, ⊤) →+* (𝓞 F)) :=
    fun J => Ideal.comap_symm e
  rw [hcomap, hcomap, hcomap, comap_mul, key, key, key, map_mul, Ideal.map_mul]

/-! ## ★★★`Proposition 1.4, (i)` -/

/-- ★★★**[GenEll] Proposition 1.4, (i)** —— 構成した高さの加法性(**無条件**)。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

★★★**仮説は原文自身のもの(`x_F` が因子を通らない)だけである。**

★★posit された `HeightTheoryData` の上ではこの主張は**偽**だった
(`Check/GenEll/HeightAxiomGap.lean` に反例)。
★★★**構成に置き換えると真になる**——それが本ファイルの到達点である。 -/
theorem htArith_tensor_unconditional {X : Scheme.{0}} (D E : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X)
    (hD : pullbackIdeal F D.divisor xF ≠ 0) (hE : pullbackIdeal F E.divisor xF ≠ 0) :
    htArith F (D.tensor E) xF = htArith F D xF + htArith F E xF :=
  htArith_tensor F D E xF (pullbackMul_all F xF) hD hE

/-- ★★**引き戻しそのものが加法的**(無条件)。

★`ADiv F` の水準での加法性——高さはこれの次数である。 -/
theorem pullbackADiv_tensor_unconditional {X : Scheme.{0}} (D E : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X)
    (hD : pullbackIdeal F D.divisor xF ≠ 0) (hE : pullbackIdeal F E.divisor xF ≠ 0) :
    pullbackADiv F (D.tensor E) xF = pullbackADiv F D xF + pullbackADiv F E xF :=
  pullbackADiv_tensor F D E xF (pullbackMul_all F xF) hD hE

/-! ## ★出典の紐付け(`.src`)

★条つきである。`Proposition 1.4` 全体には (i)〜(iv) の 4 条すべてが要り、
(iv)(Northcott)は未着手である。 -/

def htArith_tensor_unconditional.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (i)(因子表示の範囲——x が D を通らない場合。原文は可逆層で X(ℚ̄) 全体)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
