import ABC3.Found.GenEll.GcdDecomp

/-!
# GenEll 第 338 ブロック —— **★★★★★★gcd 分解は全単射である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★和の入れ替えに必要な形にした

第 337 で「零でない格子点は `d·w`(`d > 0`、`w` 原始)に一意に分解する」を取った。
★本ブロックはそれを **`Equiv`**(全単射)に仕立てる:

    DecompIndex := {(d, w) | d > 0, gcd w = 1}   ≃   {v ∈ ℤ² | v ≠ 0}

★★これで `Σ_{v ≠ 0} f(v) = Σ_{(d,w)} f(d·w)` が `Equiv.tsum_eq` で書き換えられる。
★★★被加数の側も `summand_gcdEquiv` で `d⁻ᵏ·(w₁ + w₂τ)⁻ᵏ` に分離する。

## ★★★★★★★これで (i) に残るのは「解析」だけ

    Σ_{v ≠ 0} (v₁ + v₂τ)⁻ᵏ
      = Σ_{(d,w)} d⁻ᵏ·(w₁ + w₂τ)⁻ᵏ     ★本ブロックまでで代数は済んだ
      = (Σ_d d⁻ᵏ)·(Σ_w (w₁+w₂τ)⁻ᵏ)    ← Fubini(絶対収束)★残る
      = ζ(k)·2·E_k(τ)                   ← mathlib の正規化との突き合わせ★残る

★`E_k` の定義は `(1/2)·Σ_{gcd(v₀,v₁)=1} (v₀z + v₁)⁻ᵏ` であることを実測した
(`eisensteinSeries a k z = ∑' x : gammaSet N 1 a, eisSummand k x z`、
`gammaSet N 1 a = {v | v ≡ a ∧ gcd (v 0) (v 1) = 1}`、`eisSummand k v z = (v₀z + v₁)^(-k)`)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `DecompIndex` | ★分解の添字 |
| `decomp_ne_zero` | ★★`d·w ≠ 0` |
| `gcdEquiv` | ★★★★★★**gcd 分解は全単射** |
| `summand_gcdEquiv` | ★★★被加数の分離(全単射に沿った形) |
-/

namespace ABC3.Found.GenEll

/-! ## ★分解の添字 -/

/-- ★**分解の添字**——`(d, w)` で `d > 0`、`w` は原始。 -/
def DecompIndex : Type := {p : ℕ × (ℤ × ℤ) // 0 < p.1 ∧ Int.gcd p.2.1 p.2.2 = 1}

/-- ★★`d > 0` と `w` 原始なら `d·w ≠ 0`。 -/
theorem decomp_ne_zero (p : DecompIndex) :
    (((p.1.1 : ℤ) * p.1.2.1, (p.1.1 : ℤ) * p.1.2.2) : ℤ × ℤ) ≠ (0, 0) := by
  intro h
  have hd : (p.1.1 : ℤ) ≠ 0 := by exact_mod_cast p.2.1.ne'
  have h1 : (p.1.1 : ℤ) * p.1.2.1 = 0 := congrArg Prod.fst h
  have h2 : (p.1.1 : ℤ) * p.1.2.2 = 0 := congrArg Prod.snd h
  have hw1 : p.1.2.1 = 0 := (mul_eq_zero.1 h1).resolve_left hd
  have hw2 : p.1.2.2 = 0 := (mul_eq_zero.1 h2).resolve_left hd
  have h3 := p.2.2
  rw [hw1, hw2] at h3
  simp at h3

/-! ## ★★★★★★全単射 -/

/-- ★★★★★★**gcd 分解は全単射である**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★単射性は第 337 の一意性、全射性は存在から出る。 -/
noncomputable def gcdEquiv : DecompIndex ≃ {v : ℤ × ℤ // v ≠ (0, 0)} :=
  Equiv.ofBijective
    (fun p => ⟨((p.1.1 : ℤ) * p.1.2.1, (p.1.1 : ℤ) * p.1.2.2), decomp_ne_zero p⟩)
    ⟨by
      rintro ⟨p, hp⟩ ⟨q, hq⟩ hpq
      have hv : (((p.1 : ℤ) * p.2.1, (p.1 : ℤ) * p.2.2) : ℤ × ℤ) ≠ (0, 0) :=
        decomp_ne_zero ⟨p, hp⟩
      obtain ⟨r, _, hru⟩ := exists_unique_gcd_decomp _ hv
      have h1 : p = r := hru p ⟨hp.1, hp.2, rfl⟩
      have h2 : q = r := by
        refine hru q ⟨hq.1, hq.2, ?_⟩
        exact congrArg Subtype.val hpq
      exact Subtype.ext (h1.trans h2.symm),
     by
      rintro ⟨v, hv⟩
      obtain ⟨r, ⟨hr1, hr2, hr3⟩, _⟩ := exists_unique_gcd_decomp v hv
      exact ⟨⟨r, hr1, hr2⟩, Subtype.ext hr3.symm⟩⟩

/-- ★★★**被加数の分離**(全単射に沿った形)。 -/
theorem summand_gcdEquiv (t : ℂ) (k : ℕ) (p : DecompIndex) :
    ((((gcdEquiv p).1.1 : ℂ) + ((gcdEquiv p).1.2 : ℂ) * t) ^ k)⁻¹
      = ((p.1.1 : ℂ) ^ k)⁻¹ * (((p.1.2.1 : ℂ) + (p.1.2.2 : ℂ) * t) ^ k)⁻¹ :=
  eisSummand_gcd_factor t p.1.1 p.1.2 k

/-! ## ★出典の紐付け(`.src`) -/

def gcdEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def summand_gcdEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
