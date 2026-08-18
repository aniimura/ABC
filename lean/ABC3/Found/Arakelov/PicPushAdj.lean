import ABC3.Found.Arakelov.PicPullIdealHom

/-!
# Arakelov (B2) 第 208 ブロック —— **引き戻した切断は比較射の像に入る**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★§9-225 の補題 B の第 1 歩

    f^# s  =  pullIdealHom.app (f⁻¹V) (unit.app V s)

★★これが「比較射の像が生成元を含む」ことの根拠である。
`f^*` の**具体形を一度も使わない**——随伴の関係式

    homEquiv g = unit ≫ pushforward.map g

だけで出るのが要点である(§9-221 で「具体形は塞がっている」と測った迂回)。

## ★★`Adjunction.homEquiv_apply` は定義的

`pullIdealHom := homEquiv.symm (pushHom)` と定義したので、
`homEquiv (homEquiv.symm x) = x`(`Equiv.apply_symm_apply`)を当てるだけで
`pushHom = unit ≫ pushforward.map pullIdealHom` が出る。
★mathlib の `homEquiv_apply` が `unit ≫ map` に**定義的に等しい**ためである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `pushHom_eq` | ★★★随伴の関係式(射のレベル) |
| `pushHom_app` | ★★★★**切断レベル——引き戻した切断は像に入る** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y) (D : Y.IdealSheafData)

/-- ★随伴の関係式:押し出し側の射は単位と比較射の合成である。 -/
theorem pushHom_eq :
    pushHom f D = (PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.app
      (idealPresheaf D)
      ≫ (PresheafOfModules.pushforward (pullbackPhi f)).map (pullIdealHom f D) := by
  rw [pullIdealHom]
  exact ((PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).homEquiv
    (idealPresheaf D) (idealPresheaf (D.comap f))).apply_symm_apply (pushHom f D) |>.symm




/-- ★★切断レベル:引き戻した切断は比較射の像に入る。 -/
theorem pushHom_app (V : Y.Opens) (s : ((idealPresheaf D).obj (op V) : Type u)) :
    ((pushHom f D).app (op V)).hom s
      = ((pullIdealHom f D).app (op (f ⁻¹ᵁ V))).hom
        ((((PresheafOfModules.pullbackPushforwardAdjunction
          (pullbackPhi f)).unit.app (idealPresheaf D)).app (op V)).hom s) :=
  congrArg (fun (m : idealPresheaf D ⟶ (PresheafOfModules.pushforward (pullbackPhi f)).obj
      (idealPresheaf (D.comap f))) =>
    (ModuleCat.Hom.hom (PresheafOfModules.Hom.app m (op V))) s) (pushHom_eq f D)


/-! ## ★出典の紐付け(`.src`) -/

def pushHom_app.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——引き戻した切断は比較射の像に入る)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
