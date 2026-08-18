import ABC3.Found.Arakelov.PicPullImage

/-!
# Arakelov (B2) 第 212 ブロック —— **像の閉性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★span 帰納に要る 3 つの閉性

`Submodule.span_induction` で「像 ⊇ 生成元 ⟹ 像 ⊇ イデアル」を回すには
**零・加法・スカラー倍**の 3 つが要る。第 209 でスカラー倍は取ったので、残り 2 つを取る。

★どれも「線型写像の像」であることから直に出る。

## ★★これで span 帰納の**枠が閉じた**

    Submodule.span_induction
      生成元  ← ★第 208(引き戻した切断は像に入る)
      零      ← ★本ブロック
      加法    ← ★本ブロック
      スカラー ← ★第 209

★★実測(2026-08-19)で、この 4 つを与えれば `span_induction` が**通る**ことを確認した。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `imgIdeal_zero` | ★零は像に入る |
| `imgIdeal_add` | ★★像は加法で閉じる |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y) (D : Y.IdealSheafData)

/-- ★**零は像に入る**。 -/
theorem imgIdeal_zero (A : X.Opens) : (0 : (Γ(X, A) : Type u)) ∈ imgIdeal f D A :=
  ⟨_, ⟨0, map_zero _⟩, rfl⟩

/-- ★★**像は加法で閉じる**。 -/
theorem imgIdeal_add (A : X.Opens) {a b : (Γ(X, A) : Type u)}
    (ha : a ∈ imgIdeal f D A) (hb : b ∈ imgIdeal f D A) : a + b ∈ imgIdeal f D A := by
  obtain ⟨x, ⟨y, rfl⟩, rfl⟩ := ha
  obtain ⟨x', ⟨y', rfl⟩, rfl⟩ := hb
  exact ⟨_, ⟨y + y', map_add _ _ _⟩, rfl⟩

/-- ★**同型で写した所属**。 -/
theorem mem_comap_iff {R S : CommRingCat.{u}} (e : R ≅ S) (I : Ideal (R : Type u))
    (x : (R : Type u)) : e.hom.hom x ∈ I.comap (e.inv.hom) ↔ x ∈ I := by
  have hx : e.inv.hom (e.hom.hom x) = x :=
    congrArg (fun (m : _ ⟶ _) => (CommRingCat.Hom.hom m) x) e.hom_inv_id
  simp only [Ideal.mem_comap, hx]

/-! ## ★出典の紐付け(`.src`) -/

def imgIdeal_add.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——像の閉性)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
