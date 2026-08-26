import ABC3.Found.GenEll.LatticeNorm
import Mathlib.AlgebraicGeometry.EllipticCurve.VariableChange

/-!
# GenEll 第 335 ブロック —— **★★★★★★束の曲線と `j` 不変量**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★束の側と曲線の側を繋ぐ

`℘'² = 4℘³ − g₂℘ − g₃` で `x = ℘`、`y = ℘'/2` と置くと `y² = x³ − (g₂/4)x − (g₃/4)`。
★これを `latticeCurve P` と呼ぶ。

★★本ブロックはその**判別式と `j` 不変量**を計算する:

    Δ(latticeCurve P) = g₂³ − 27g₃² = latticeDisc P
    j(latticeCurve P) = 1728·g₂³ / (g₂³ − 27g₃²)

★★★どちらも古典的な式そのものであり、`b₂ b₄ b₆ b₈` を展開して `ring` で出る。

## ★★★★これで (i) の位置づけが確定する

`latticeCurve P` が**楕円曲線である**(`Δ` が単元)ことは
**`latticeDisc P ≠ 0`** とちょうど同値である(`isElliptic_latticeCurve` と逆向き)。
★したがって一意化の連鎖の (i) は「`latticeCurve` が楕円曲線であること」に他ならない。

★★`j` の式は (iii)「`j` の全射性」を格子の言葉に翻訳する:
与えられた `j₀` に対し `1728g₂³/(g₂³−27g₃²) = j₀` なる束を作ればよい。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `latticeCurve` | ★★束に対応する Weierstrass 曲線 |
| `latticeCurve_Δ` | ★★★★★**`Δ = g₂³ − 27g₃²`** |
| `isElliptic_latticeCurve` | ★★★★判別式が非零なら楕円曲線 |
| `latticeDisc_ne_zero_of_isElliptic` | ★★★逆向き |
| `latticeCurve_c₄` | ★★`c₄ = 12g₂` |
| `latticeCurve_j` | ★★★★★★**`j = 1728g₂³/(g₂³−27g₃²)`** |
-/

namespace ABC3.Found.GenEll

open PeriodPair WeierstrassCurve

/-! ## ★★束に対応する Weierstrass 曲線 -/

/-- ★★**束 `Λ` に対応する Weierstrass 曲線** `y² = x³ − (g₂/4)x − (g₃/4)`。

★`℘'² = 4℘³ − g₂℘ − g₃` で `x = ℘`、`y = ℘'/2` と置いた形である。 -/
noncomputable def latticeCurve (P : PeriodPair) : WeierstrassCurve ℂ :=
  ⟨0, 0, 0, -P.g₂ / 4, -P.g₃ / 4⟩

/-- ★★★★★**判別式は `g₂³ − 27g₃²` である**。 -/
theorem latticeCurve_Δ (P : PeriodPair) : (latticeCurve P).Δ = latticeDisc P := by
  simp only [latticeCurve, latticeDisc, WeierstrassCurve.Δ, WeierstrassCurve.b₂,
    WeierstrassCurve.b₄, WeierstrassCurve.b₆, WeierstrassCurve.b₈]
  ring

/-- ★★★★判別式が非零なら `latticeCurve` は楕円曲線である。 -/
theorem isElliptic_latticeCurve (P : PeriodPair) (h : latticeDisc P ≠ 0) :
    (latticeCurve P).IsElliptic :=
  ⟨by rw [latticeCurve_Δ]; exact isUnit_iff_ne_zero.2 h⟩

/-- ★★★逆向き——楕円曲線なら判別式は非零。 -/
theorem latticeDisc_ne_zero_of_isElliptic (P : PeriodPair) [h : (latticeCurve P).IsElliptic] :
    latticeDisc P ≠ 0 := by
  have h2 := h.isUnit
  rw [latticeCurve_Δ] at h2
  exact h2.ne_zero

/-! ## ★★★★★★`j` 不変量 -/

/-- ★★`c₄ = 12g₂`。 -/
theorem latticeCurve_c₄ (P : PeriodPair) : (latticeCurve P).c₄ = 12 * P.g₂ := by
  simp only [latticeCurve, WeierstrassCurve.c₄, WeierstrassCurve.b₂, WeierstrassCurve.b₄]
  ring

/-- ★★★★★★**`j = 1728·g₂³/(g₂³ − 27g₃²)`**——古典的な式。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★これが (iii)「`j` の全射性」を格子の言葉に翻訳する。 -/
theorem latticeCurve_j (P : PeriodPair) [(latticeCurve P).IsElliptic] :
    (latticeCurve P).j = 1728 * P.g₂ ^ 3 / latticeDisc P := by
  rw [WeierstrassCurve.j, latticeCurve_c₄, Units.val_inv_eq_inv_val]
  have hD : ((latticeCurve P).Δ' : ℂ) = latticeDisc P := by
    rw [WeierstrassCurve.coe_Δ', latticeCurve_Δ]
  rw [hD]
  have hne : latticeDisc P ≠ 0 := latticeDisc_ne_zero_of_isElliptic P
  field_simp
  ring

/-! ## ★出典の紐付け(`.src`) -/

def latticeCurve.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def latticeCurve_j.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
