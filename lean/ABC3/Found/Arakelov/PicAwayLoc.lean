import ABC3.Found.Arakelov.PicAwayTower

/-!
# Arakelov (B1) 第 96 ブロック —— **`R_{t·g}` は `R_g` の局所化**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★mathlib の向きと逆だったので回った

★§9-99 で測った通り、mathlib の在庫は**すべて逆向き**である:

| 在庫 | 向き |
|---|---|
| `IsLocalization.Away.mul` | `Away x S` ∧ `Away (algebraMap y) T` ⟹ `Away (y·x) T` |
| `isLocalization_of_submonoid_le` | `M ≤ N` を要求(`powers g ≤ powers (t·g)` は**偽**) |

★★**回り道**: `powers (t·g)` の代わりに **`closure {t, g}`** を使う。
これなら `powers g ≤ closure {t,g}` が成り立つ。

## ★★機構 —— 2 手

| 手 | 内容 |
|---|---|
| 1 | `powers (t·g) ≤ closure {t,g}` かつ `closure` の元は `R_{t·g}` で可逆 → `IsLocalization.of_le` |
| 2 | `powers g ≤ closure {t,g}` → `isLocalization_of_submonoid_le` |

★★★**第 3 手(`powers (algebraMap t)` へ絞る)は不要**である
——`Module.free_of_isLocalizedModule` は**任意の積閉集合**で使えるので、
`closure {t,g}` の像のままで良い(2026-08-18 実測)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isLocalization_closure` | ★`R_{t·g}` は `closure {t,g}` での局所化 |
| `isLocalization_closure_map` | ★★★★**`R_{t·g}` は `R_g` の局所化** |
-/

universe u

namespace ABC3.Found.Arakelov

variable (R : Type u) [CommRing R] (g t : R)

/-- ★**`R_{t·g}` は `closure {t, g}` での局所化である**。 -/
theorem isLocalization_closure : IsLocalization (Submonoid.closure ({t, g} : Set R))
    (Localization (Submonoid.powers (t * g))) := by
  refine IsLocalization.of_le (Submonoid.powers (t * g)) _ ?_ ?_
  · rw [Submonoid.powers_le]
    exact Submonoid.mul_mem _ (Submonoid.subset_closure (by simp))
      (Submonoid.subset_closure (by simp))
  · intro r hr
    induction hr using Submonoid.closure_induction with
    | mem x hx =>
      rcases hx with h | h <;> rw [h]
      · exact IsLocalization.Away.isUnit_of_dvd (x := t * g) ⟨g, rfl⟩
      · exact IsLocalization.Away.isUnit_of_dvd (x := t * g) ⟨t, by ring⟩
    | one => exact isUnit_one.map _
    | mul x y _ _ hx hy => rw [map_mul]; exact hx.mul hy

/-- ★★★★**`R_{t·g}` は `R_g` の局所化である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで自由性を `D(g)` から `D(t·g)` へ運べる。 -/
theorem isLocalization_closure_map :
    IsLocalization ((Submonoid.closure ({t, g} : Set R)).map
        (algebraMap R (Localization (Submonoid.powers g))))
      (Localization (Submonoid.powers (t * g))) := by
  haveI := isLocalization_closure R g t
  exact IsLocalization.isLocalization_of_submonoid_le
    (Localization (Submonoid.powers g)) (Localization (Submonoid.powers (t * g)))
    (Submonoid.powers g) (Submonoid.closure ({t, g} : Set R))
    (by rw [Submonoid.powers_le]; exact Submonoid.subset_closure (by simp))

/-! ## ★出典の紐付け(`.src`) -/

def isLocalization_closure_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——R_{t·g} は R_g の局所化)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
