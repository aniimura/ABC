import ABC3.Found.Arakelov.ArcAppIsoNat

/-!
# Arakelov (C3) 第 268 ブロック —— ★★★★★★**橋の自然性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★これで橋の 3 法則が揃った

| 法則 | ブロック |
|---|---|
| 加法性 | 第 255 |
| スカラー両立 | ★第 266 |
| ★**自然性** | ★★**本ブロック** |

★★これで `bridgeApp` は前層加群の射になり、`bridgeApp_inv`(第 255)と合わせて
**`restrict F V.ι ≅ 𝒪_V`(§9-297 の橋)**が組める。

## ★★機構 —— `rfl` 補題 2 本 + `e.hom` の自然性 + `appIso` の自然性

| 段 | 使うもの |
|---|---|
| 合成の適用 | ★`msplit`(**明示束縛子の `rfl`**、4 例目) |
| 単位対象の制限射 | ★`unitMap_val`(`rfl`) |
| `e.hom` の自然性 | ★mathlib(`Hom.naturality`) |
| `appIso` の自然性 | ★第 267 |

★★★`e.hom.naturality` の片側には **`restrictScalars` が噛んでいる**が、
元素レベルでは同じ関数なので `msplit` で割れば消える。

## ★摩擦(再掲)—— `rw ... at h` は二重路で落ちる

`rw [msplit, msplit] at hnat` は `presheaf` の二重路で落ちた。
★`(msplit …).symm.trans (hnat.trans (msplit …))` と**項で書くと通る**。
★★これで「`rw` をやめて `trans` で繋ぐ」は本セッション **5 例目**である。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{0}} (V : X.Opens) (F : X.Modules)

/-- ★`ModuleCat` の合成の適用。 -/
theorem msplit {R : RingCat.{0}} (A B C : ModuleCat.{0} R) (f : A ⟶ B) (g : B ⟶ C)
    (x : (A : Type)) : ((f ≫ g).hom) x = g.hom (f.hom x) := rfl

/-- ★★単位対象の制限射は関数の制限である。 -/
theorem unitMap_val {W W' : V.toScheme.Opens} (k : W' ⟶ W)
    (y : (((𝟙_ (PresheafModulesOn X V)).obj (op (overObj V W))) : Type)) :
    unitValAt V W' (((𝟙_ (PresheafModulesOn X V)).map (overMap V k).op).hom y)
      = (X.presheaf.map ((Scheme.Hom.opensFunctor V.ι).map k).op).hom (unitValAt V W y) :=
  rfl


variable (e : (restrictPresheafFunctor X V).obj F.val ≅ 𝟙_ (PresheafModulesOn X V))

/-- ★★★★橋の自然性。 -/
theorem bridgeApp_nat {W W' : (V.toScheme.Opens)ᵒᵖ} (i : W ⟶ W')
    (x : ((Scheme.Modules.restrict F V.ι).val.obj W : Type)) :
    bridgeApp F V e W' (((Scheme.Modules.restrict F V.ι).val.map i).hom x)
      = ((unitModules V.toScheme).val.map i).hom (bridgeApp F V e W x) := by
  have hnat := congrArg (fun (m : _ ⟶ _) => (ModuleCat.Hom.hom m) (secOverV V F W.unop x))
    (e.hom.naturality (overMap V i.unop).op)
  have hnat' : (e.hom.app (op (overObj V W'.unop))).hom
      ((((restrictPresheafFunctor X V).obj F.val).map (overMap V i.unop).op).hom
        (secOverV V F W.unop x))
      = ((𝟙_ (PresheafModulesOn X V)).map (overMap V i.unop).op).hom
        ((e.hom.app (op (overObj V W.unop))).hom (secOverV V F W.unop x)) :=
    (msplit _ _ _ _ _ _).symm.trans (hnat.trans (msplit _ _ _ _ _ _))
  show (V.ι.appIso W'.unop).hom.hom (unitValAt V W'.unop
      ((e.hom.app (op (overObj V W'.unop))).hom
        ((((restrictPresheafFunctor X V).obj F.val).map (overMap V i.unop).op).hom
          (secOverV V F W.unop x)))) = _
  rw [hnat', unitMap_val]
  exact congrArg (fun (m : _ ⟶ _) => (CommRingCat.Hom.hom m)
    (unitValAt V W.unop ((e.hom.app (op (overObj V W.unop))).hom (secOverV V F W.unop x))))
    (appIso_hom_naturality V i).symm


/-! ## ★出典の紐付け(`.src`) -/

def bridgeApp_nat.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——橋の自然性)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
