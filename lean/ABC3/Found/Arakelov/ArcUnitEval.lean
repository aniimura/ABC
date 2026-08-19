import ABC3.Found.Arakelov.ArcEvalNat

/-!
# Arakelov (C3) 第 251 ブロック —— ★★★★★**構造層のファイバー評価は関数の値**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★連続性の最後の橋

自明化 `t : L|_V ≅ 𝒪_V` を使うと、ノルムは

    ‖s‖(p) = |(t s)(p)|

になる。★右辺が**正則関数の点での値**であることを言うのが本ブロックである:

    (unitFiberIso p).hom (arcEval p 𝒪_X f) = p^♯(f)        (`evalUnit_eq`)

★★これで `continuous_evalAffine`(第 5 ブロック)に繋がる。

## ★★機構 —— mathlib の随伴の特徴づけを `⊤` に落とすだけ

| 段 | 使うもの |
|---|---|
| `pullbackObjUnitToUnit` の特徴づけ | ★mathlib `pullbackPushforwardAdjunction_homEquiv_pullbackObjUnitToUnit` |
| `homEquiv f = η ≫ G.map f` | ★mathlib `Adjunction.homEquiv_unit` |
| `unitToPushforwardObjUnit` の切断レベルの値 | ★mathlib `unitToPushforwardObjUnit_val_app_apply` |
| 押し出しは `⊤` では何もしない | ★`rfl`(`p ⁻¹ᵁ ⊤ = ⊤`) |

★★★**3 つとも mathlib に在った**——`Sheaf/PullbackFree.lean` である。
★在庫を測ってから書いたので、証明は 4 行で済んだ。

## ★摩擦 —— `rw [hsplit]` は要らなかった

合成の分解は `rfl` なので、`rw` すると「パターンが見つからない」で落ちる。
★★**`rfl` な補題は `rw` しない**——`exact` に任せる。
これは第 250 の教訓(`Eq.trans` は reducible しか使わない)の裏返しである。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite

variable {X : Scheme.{0}}
noncomputable def unitFiberIso (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    arcFiber p (unitModules X)
      ≅ moduleSpecΓFunctor (R := CommRingCat.of ℂ) |>.obj
          (unitModules (Spec (CommRingCat.of ℂ))) :=
  (moduleSpecΓFunctor (R := CommRingCat.of ℂ)).mapIso (pullbackUnitIso p)

/-- ★随伴の特徴づけを ⊤ の切断に落とす。 -/
theorem adj_at_top (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (f : ((X.presheaf.obj (op ⊤)) : Type)) :
    (((Scheme.Modules.pullbackPushforwardAdjunction p).homEquiv _ _
        (SheafOfModules.pullbackObjUnitToUnit p.toRingCatSheafHom)).val.app (op ⊤)).hom f
      = ((SheafOfModules.unitToPushforwardObjUnit p.toRingCatSheafHom).val.app (op ⊤)).hom f :=
  congrArg (fun (m : unitModules X ⟶ (Scheme.Modules.pushforward p).obj
      (unitModules (Spec (CommRingCat.of ℂ)))) => (m.val.app (op ⊤)).hom f)
    (SheafOfModules.pullbackPushforwardAdjunction_homEquiv_pullbackObjUnitToUnit
      p.toRingCatSheafHom)

/-- ★押し出しは `⊤` では何もしない。 -/
theorem pushforward_at_top (p : Spec (CommRingCat.of ℂ) ⟶ X)
    {N N' : (Spec (CommRingCat.of ℂ)).Modules} (ψ : N ⟶ N')
    (w : ((((Scheme.Modules.pushforward p).obj N).val.obj (op ⊤)) : Type)) :
    (((Scheme.Modules.pushforward p).map ψ).val.app (op ⊤)).hom w
      = (ψ.val.app (op ⊤)).hom w :=
  rfl

/-- ★★★构造層のファイバーでの評価は関数の値である。 -/
theorem evalUnit_eq (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (f : ((X.presheaf.obj (op ⊤)) : Type)) :
    (unitFiberIso p).hom.hom (arcEval p (unitModules X) f)
      = (Scheme.Hom.toRingCatSheafHom p).hom.app (op ⊤) f := by
  have h := adj_at_top p f
  rw [Adjunction.homEquiv_unit] at h
  rw [SheafOfModules.unitToPushforwardObjUnit_val_app_apply] at h
  exact h


/-! ## ★出典の紐付け(`.src`) -/

def evalUnit_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——構造層のファイバー評価は関数の値)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
