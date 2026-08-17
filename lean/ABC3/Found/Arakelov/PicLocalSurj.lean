import ABC3.Found.Arakelov.PicRestrictTensor

/-!
# Arakelov (B1) 第 6 ブロック —— **局所全射性はテンソルで保たれる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★結合律の残り 2 条のうち 1 条

結合律は `W (η ▷ C)`(`η` は層化の単位)に帰着し、
`WEqualsLocallyBijective` により **局所全射 + 局所単射**の 2 条になる。

★★★**局所全射の側は局所自由性を使わずに直接証明できる。**

## ★機構

`x ∈ (Q ⊗ M)(U)` は純テンソルの有限和である。
★各因子 `q` について `imageSieve f q` は被覆篩なので、
有限個の交わりも被覆篩(`J.intersection_covering`)。

★★証明の型は **`TensorProduct.induction_on`**:

| 場合 | 使う道具 |
|---|---|
| `0` | 篩は `⊤`(`0` が原像) |
| `q ⊗ₜ m` | `Presheaf.imageSieve_mem`(`f` の局所全射性) |
| `s + t` | `J.intersection_covering` |

★いずれも `J.superset_covering` で目的の篩へ落とす。
-/

universe u

namespace ABC3.Found.Arakelov

open CategoryTheory MonoidalCategory Opposite

/-- ★★★**局所全射な射をテンソルしても局所全射である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが結合律の 2 条のうち 1 条である。**局所自由性は要らない。**

★機構は `TensorProduct.induction_on` と `J.intersection_covering`。 -/
theorem isLocallySurjective_whiskerRight {C : Type u} [Category.{u} C]
    (J : GrothendieckTopology C) {R : Cᵒᵖ ⥤ CommRingCat.{u}}
    {P Q : PresheafOfModules.{u} (R ⋙ forget₂ CommRingCat.{u} RingCat.{u})}
    (f : P ⟶ Q) [PresheafOfModules.IsLocallySurjective J f]
    (M : PresheafOfModules.{u} (R ⋙ forget₂ CommRingCat.{u} RingCat.{u})) :
    PresheafOfModules.IsLocallySurjective J (f ▷ M) := by
  constructor
  intro U s
  show Presheaf.imageSieve _ s ∈ J U
  induction s using TensorProduct.induction_on with
  | zero =>
      refine J.superset_covering ?_ (J.top_mem U)
      intro V i _
      exact ⟨0, (map_zero _).trans (map_zero _).symm⟩
  | tmul q m =>
      refine J.superset_covering ?_
        (Presheaf.imageSieve_mem J ((PresheafOfModules.toPresheaf _).map f) q)
      rintro V i ⟨p, hp⟩
      refine ⟨p ⊗ₜ M.map i.op m, ?_⟩
      have hp' : f.app (op V) p = Q.map i.op q := hp
      change f.app (op V) p ⊗ₜ M.map i.op m = Q.map i.op q ⊗ₜ M.map i.op m
      rw [hp']
      rfl
  | add s t hs ht =>
      refine J.superset_covering ?_ (J.intersection_covering hs ht)
      rintro V i ⟨⟨a, ha⟩, ⟨b, hb⟩⟩
      refine ⟨a + b, ?_⟩
      rw [map_add, ha, hb, ← map_add]
      rfl

/-! ## ★出典の紐付け(`.src`) -/

def isLocallySurjective_whiskerRight.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——局所全射性がテンソルで保たれること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
