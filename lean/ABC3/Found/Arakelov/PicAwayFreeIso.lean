import ABC3.Found.Arakelov.PicAwayFreeProp

/-!
# Arakelov (B1) 第 100 ブロック —— **基底の全体で `M_h ≅ R_h`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★局所自明性の切断側が揃った

第 76 ブロック: 各素点の近傍に `M_r ≅ R_r` となる `D(r)` がある。
第 99 ブロック: 自由性は `D(g)` から `D(t·g)` へ運べる。

★★合わせると、**`D(g)` の基底 `{D(t·g)}` の全体で `M_{t·g} ≅ R_{t·g}`** が成り立つ。

## ★★機構

★`Module.Free` + `Module.Invertible` ⟹ `≃ₗ` は
mathlib の `Module.Invertible.free_iff_linearEquiv`。
★★`M_{t·g}` が可逆であることは mathlib の instance
(`Module.Invertible (Localization S) (LocalizedModule S M)`)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `linearEquiv_away_mul` | ★★★★**基底の全体で `M_{t·g} ≅ R_{t·g}`** |

## ★★★残る構造上の 1 点(記録)

★`IsLocallyTrivial` は `(restrict V).obj P ≅ 𝟙_`(**`Over V` 上の前層加群の同型**)を要求する。
★★切断ごとの同型は本ブロックで**基底の上では**出たが、
`Over V` の対象は `V` 以下の**すべての開集合**である。

★★★**「基底で同型 ⟹ 全体で同型」を `Over V` の site で言う器具が要る**
——第 90 ブロックは `TopCat.Presheaf`(空間の上)の版であり、そのままでは当たらない。

★候補: (a) 第 90 を `Over V` へ移す (b) 茎で言う(第 77 の `Over V` 版)。
-/

universe u

namespace ABC3.Found.Arakelov

variable (R : Type u) [CommRing R] (g t : R)
  (M : Type u) [AddCommGroup M] [Module R M] [Module.Invertible R M]

/-- ★★★★**基底の全体で `M_{t·g} ≅ R_{t·g}`**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★第 76(近傍で自由)と第 99(自由性の伝播)の合成である。 -/
theorem linearEquiv_away_mul
    [Module.Free (Localization (Submonoid.powers g))
      (LocalizedModule (Submonoid.powers g) M)] :
    Nonempty (LocalizedModule (Submonoid.powers (t * g)) M
      ≃ₗ[Localization (Submonoid.powers (t * g))] Localization (Submonoid.powers (t * g))) := by
  haveI := free_away_mul R g t M
  exact Module.Invertible.free_iff_linearEquiv.mp inferInstance

/-! ## ★出典の紐付け(`.src`) -/

def linearEquiv_away_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——基底の全体で M_h ≅ R_h)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
