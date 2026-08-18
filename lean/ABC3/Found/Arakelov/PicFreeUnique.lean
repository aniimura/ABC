import ABC3.Found.Arakelov.PicUnitApp

/-!
# Arakelov (B1) 第 105 ブロック —— **一点添字の自由加群は階数 1**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★`iso.inv` を計算せずに済ませる

§9-114 で

    (unitHomOfSection s).app W = (iso.inv.app W) ≫ ((freeYonedaEquiv.symm s).app W)

の **`iso.inv` の値を計算する必要は無い**(同型だから全単射)と気づいた。
★★残るのは右側——`free(Hom(W,T)) → P(W)` の全単射性である。

★★★`Hom(W,T)` は**一点**(終対象)なので、`free` は**階数 1**である。

## ★★詰まりの回避 —— `ModuleCat.freeMk` を `Finsupp.single` に書き直す

★`ModuleCat.freeMk x` のままだと `LinearEquiv.map_smul` が**当たらない**
(`ModuleCat` の `•` と `Finsupp` の `•` が別経路)。
★★**`Finsupp.single x 1` に書き直すと `simp` が通る**(2026-08-24 実測)。
★★★両者は **`rfl`** である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `freeMk_eq_single` | ★`freeMk x = Finsupp.single x 1`(`rfl`) |
| `uniqueLinearEquiv_symm_apply` | ★★★**階数 1 の同型の逆は「生成元の倍」** |
-/

universe u

namespace ABC3.Found.Arakelov

open CategoryTheory

variable {A : Type u} [CommRing A] {ι : Type u} [Unique ι]

/-- ★**`freeMk x` は `Finsupp.single x 1` である**(`rfl`)。

★★これに書き直すと `LinearEquiv` の補題が当たるようになる。 -/
theorem freeMk_eq_single (x : ι) :
    (ModuleCat.freeMk x : (ModuleCat.free A).obj ι) = Finsupp.single x 1 := rfl

/-- ★★★**階数 1 の同型の逆は「生成元の倍」である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで `free(一点) → P(W)` の全単射性が
「生成元の像による乗法」の全単射性(第 103 ブロック)に落ちる。 -/
theorem uniqueLinearEquiv_symm_apply (c : A) :
    (Finsupp.uniqueLinearEquiv A A (default : ι)).symm c
      = c • (Finsupp.single (default : ι) (1 : A)) := by
  apply (Finsupp.uniqueLinearEquiv A A (default : ι)).injective
  simp [Finsupp.uniqueLinearEquiv]

/-! ## ★出典の紐付け(`.src`) -/

def uniqueLinearEquiv_symm_apply.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——一点添字の自由加群は階数 1)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
