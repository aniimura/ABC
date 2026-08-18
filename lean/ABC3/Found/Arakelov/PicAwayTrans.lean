import ABC3.Found.Arakelov.PicUnitHom

/-!
# Arakelov (B1) 第 94 ブロック —— **局所化の推移**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★局所自明性に要る「小さくする」段

第 76 ブロックで「可逆加群は `D(r)` の上で自由」が出た。
★★局所自明性を示すには `D(r)` の**部分**でも自由でなければならない。

★★★`D(t·r) = D(t) ⊓ D(r)` は `D(r)` の**基底**をなすので、

    R_{t·r} は R_r の(t の像での)局所化である

を示せば、`Module.free_of_isLocalizedModule` で `M_{t·r}` の自由性が出る。

## ★★機構 —— mathlib の在庫 2 つ

| 在庫 | 内容 |
|---|---|
| `IsLocalization.Away.isUnit_of_dvd` | ★`g ∣ t·g` なので `g` は `R_{t·g}` で可逆 |
| `IsLocalization.Away.lift` | ★環準同型 `R_g → R_{t·g}` |
| `IsLocalization.Away.mul` | ★★`R_{t·g}` は `R_g` の `t` での局所化 |

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `awayToAwayMul` | ★★**環準同型 `R_g → R_{t·g}`** |
-/

universe u

namespace ABC3.Found.Arakelov

variable (R : Type u) [CommRing R] (g t : R)

/-- ★★**環準同型 `R_g → R_{t·g}`**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★`g ∣ t·g` なので `g` は `R_{t·g}` の中で可逆であり、普遍性で持ち上がる。 -/
noncomputable def awayToAwayMul : Localization (Submonoid.powers g) →+*
    Localization (Submonoid.powers (t * g)) :=
  IsLocalization.Away.lift (S := Localization (Submonoid.powers g)) g
    (IsLocalization.Away.isUnit_of_dvd (S := Localization (Submonoid.powers (t * g)))
      (x := t * g) ⟨t, by ring⟩)

/-! ## ★出典の紐付け(`.src`) -/

def awayToAwayMul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——局所化の推移)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
