import ABC3.Found.GenEll.ComapLocal
import ABC3.Found.GenEll.StalkSupport

/-!
# [GenEll] Proposition 1.4, (i) —— **茎の引き戻しと `comap` を繋ぐ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

## ★★★本日作った 2 つの引き戻しを繋ぐ

本日、イデアル層の引き戻しを **2 通り**作った:

| 構成 | 場所 | 性質 |
|---|---|---|
| `stalkPullback` —— 茎ごとの押し出し | `StalkPullback.lean` | ★**積を保つ**(`stalkPullback_mul`) |
| `comap` —— mathlib のファイバー積の核 | mathlib | ★積を保つかは未証明 |

★★**この 2 つが一致すれば `comap_mul` が出る。**
本ファイルは**片側の包含**を取る:

    `stalkPullback I f x ≤ (I.comap f).stalk x`

★逆向きにはファイバー積の `Γ` の計算が要る(`ComapLocal.lean` の記録を参照)。

## ★★機構

`ComapLocal.lean` の `map_appLE_le_ideal_comap`(切断の水準での包含)を
**茎へ持ち上げる**だけである。要るのは 2 本の射の等式:

- `germ_stalkMap` —— 芽と茎写像の両立
- `germ_res` —— 芽と制限の両立

★★★**§9-13 の規則どおり、射のレベルで合成を作ってから 1 度だけ降りる。**
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory CategoryTheory.Limits

variable {X Y : Scheme}

/-! ## ★★芽・茎写像・制限の三者が繋がる -/

/-- ★★**`germ ≫ stalkMap = appLE ≫ germ`**。

★`germ_stalkMap`(芽と茎写像)と `germ_res`(芽と制限)を繋いだもの。
★★これが「切断の水準の包含」を「茎の水準」へ持ち上げる橋である。 -/
theorem germ_comp_stalkMap (f : X ⟶ Y) (x : X) (U : X.affineOpens) (hxU : x ∈ U.1)
    (V : Y.affineOpens) (h : U.1 ≤ f ⁻¹ᵁ V.1) :
    Y.presheaf.germ V.1 (f.base x) (h hxU) ≫ f.stalkMap x
      = f.appLE V.1 U.1 h ≫ X.presheaf.germ U.1 x hxU := by
  rw [Scheme.Hom.germ_stalkMap, Scheme.Hom.appLE, Category.assoc]
  congr 1
  erw [X.presheaf.germ_res]

/-! ## ★★★包含 -/

/-- ★★★**茎ごとの引き戻しは `comap` の茎に含まれる**。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

★★★**本日作った 2 つの引き戻しが繋がった。**
`stalkPullback` は積を保つ(`stalkPullback_mul`)ので、
**逆向きの包含が入れば `comap_mul` が出る。**

★機構は `map_appLE_le_ideal_comap` を `Ideal.map_mono` で茎へ持ち上げるだけ。 -/
theorem stalkPullback_le_comap_stalk (I : Y.IdealSheafData) (f : X ⟶ Y) (x : X)
    (U : X.affineOpens) (hxU : x ∈ U.1) (V : Y.affineOpens) (h : U.1 ≤ f ⁻¹ᵁ V.1) :
    stalkPullback I f x ≤ (I.comap f).stalk x := by
  rw [Scheme.IdealSheafData.stalk_eq_map (U := U) hxU, stalkPullback,
    Scheme.IdealSheafData.stalk_eq_map (U := V) (h hxU), Ideal.map_map,
    ← CommRingCat.hom_comp, germ_comp_stalkMap f x U hxU V h, CommRingCat.hom_comp,
    ← Ideal.map_map]
  exact Ideal.map_mono (map_appLE_le_ideal_comap I f U V h)

/-! ## ★出典の紐付け(`.src`) -/

def stalkPullback_le_comap_stalk.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (i)(2 つの引き戻しの片側の包含のみ)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
