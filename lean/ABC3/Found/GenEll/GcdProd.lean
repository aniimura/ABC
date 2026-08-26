import ABC3.Found.GenEll.GcdEquiv

/-!
# GenEll 第 339 ブロック —— **★★★★★★添字を「正整数 × 原始ベクトル」の積にする**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★Fubini を当てられる形にする

第 338 で `DecompIndex ≃ {v ≠ 0}` を取った。★ただし `DecompIndex` は
`{p : ℕ × (ℤ×ℤ) // 0 < p.1 ∧ gcd p.2 = 1}` という**条件の連言つき部分型**であって、
`tsum` の Fubini(`∑' (a,b), f a * g b = (∑' a, f a)(∑' b, g b)`)を当てるには
**積の形**でなければならない。

★★本ブロックはそれを

    {d : ℕ // 0 < d} × {w : ℤ × ℤ // gcd w = 1}  ≃  {v : ℤ × ℤ // v ≠ 0}

に整える(`Equiv.subtypeProdEquivProd` と同じ形を手で作る)。
★★★被加数も積の形 `d⁻ᵏ · (w₁ + w₂τ)⁻ᵏ` で書ける(`summand_nonzeroProdEquiv`)。

## ★★残り((i) の解析)

★`Σ_d d⁻ᵏ` の絶対収束(`k ≥ 2`)と
`Σ_{w 原始} (w₁+w₂τ)⁻ᵏ` の絶対収束(`k ≥ 3`)。
★★mathlib は `EisensteinSeries.summable_one_div_norm_rpow`(`2 < k` で
`Summable fun x => ‖x‖ ^ (-k)`)を持つ(2026-08-26 実測)。
★★★そこから `tsum_mul_tsum_of_summable_norm` で積に分ける。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `decompProdEquiv` | ★★連言つき部分型を積に |
| `nonzeroProdEquiv` | ★★★★★★**`{d>0} × {原始} ≃ {v ≠ 0}`** |
| `summand_nonzeroProdEquiv` | ★★★被加数が積の形になる |
-/

namespace ABC3.Found.GenEll

/-! ## ★★連言つき部分型を積に -/

/-- ★★分解の添字を積の形に整える。 -/
noncomputable def decompProdEquiv :
    DecompIndex ≃ ({d : ℕ // 0 < d} × {w : ℤ × ℤ // Int.gcd w.1 w.2 = 1}) where
  toFun p := (⟨p.1.1, p.2.1⟩, ⟨p.1.2, p.2.2⟩)
  invFun q := ⟨(q.1.1, q.2.1), q.1.2, q.2.2⟩
  left_inv := by rintro ⟨⟨d, w⟩, hd, hw⟩; rfl
  right_inv := by rintro ⟨⟨d, hd⟩, ⟨w, hw⟩⟩; rfl

/-- ★★★★★★**零でない格子点は「正の整数 × 原始ベクトル」と全単射である**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any -/
noncomputable def nonzeroProdEquiv :
    ({d : ℕ // 0 < d} × {w : ℤ × ℤ // Int.gcd w.1 w.2 = 1}) ≃ {v : ℤ × ℤ // v ≠ (0, 0)} :=
  decompProdEquiv.symm.trans gcdEquiv

/-- ★★★**被加数が積の形になる**——Fubini を当てる直前の形。 -/
theorem summand_nonzeroProdEquiv (t : ℂ) (k : ℕ)
    (q : {d : ℕ // 0 < d} × {w : ℤ × ℤ // Int.gcd w.1 w.2 = 1}) :
    ((((nonzeroProdEquiv q).1.1 : ℂ) + ((nonzeroProdEquiv q).1.2 : ℂ) * t) ^ k)⁻¹
      = ((q.1.1 : ℂ) ^ k)⁻¹ * (((q.2.1.1 : ℂ) + (q.2.1.2 : ℂ) * t) ^ k)⁻¹ :=
  summand_gcdEquiv t k (decompProdEquiv.symm q)

/-! ## ★出典の紐付け(`.src`) -/

def nonzeroProdEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def summand_nonzeroProdEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
