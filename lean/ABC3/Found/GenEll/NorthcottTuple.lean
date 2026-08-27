import ABC3.Found.GenEll.NorthcottClassical

/-!
# [GenEll] Proposition 1.4, (iv) の受け皿 —— **次数有界 Northcott の組版**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

## ★★★★なぜ「組」が要るのか

`Proposition 1.4, (iv)` は `X` 上の高さについての Northcott 性である。
★原文の証明は **`X` の射影埋め込み**を経由する——`X` が `ℤ`-固有だから
`ℙⁿ` に閉埋め込みでき、高さは `ℙⁿ` の高さと `O(1)` で一致する。

★★したがって最後に要るのは **`ℙⁿ(ℚ̄)` の次数有界 Northcott** であり、
座標で書けば「**次数 `≤ d`・高さ `≤ B` の代数的数の組は有限個**」である。

★★★1 成分の場合(`X = ℙ¹`)は `NorthcottClassical.lean` で証明済みなので、
本ファイルは**それを組に持ち上げるだけ**である(`Set.Finite.pi`)。

## ★本ファイルで取れるもの

| 定理 | 内容 |
|---|---|
| `boundedAlg` | ★次数 `≤ d`・高さ `≤ B` の代数的数の集合 |
| `boundedAlg_finite` | ★**有限**(古典的 Northcott そのもの) |
| `finite_pi_boundedAlg` | ★★**組も有限** |
| `finite_of_injOn_boundedAlg` | ★★★★**そこへ単射に写る集合は有限**——射影モデル経由で使う受け皿 |
-/

namespace ABC3.Found.GenEll

open NumberField IntermediateField

/-! ## ★次数と高さで抑えた代数的数 -/

/-- ★**`ℚ` 上の次数が `d` 以下、かつ高さが `B` 以下の代数的数**の集合。

★体を固定していない——原文の `X(ℚ̄)^{≤d}` の形である。 -/
def boundedAlg (d : ℕ) (B : ℝ) : Set ℂ :=
  {α : ℂ | ∃ h : IsIntegral ℚ α,
    Module.finrank ℚ ℚ⟮α⟯ ≤ d ∧
    (haveI := IntermediateField.adjoin.finiteDimensional h
     haveI : NumberField ℚ⟮α⟯ := ⟨⟩
     Height.mulHeight₁ (⟨α, IntermediateField.mem_adjoin_simple_self ℚ α⟩ : ℚ⟮α⟯) ≤ B)}

/-- ★**古典的 Northcott** —— `boundedAlg` は有限集合である。 -/
theorem boundedAlg_finite (d : ℕ) (B : ℝ) : (boundedAlg d B).Finite :=
  finite_of_finrank_le_of_mulHeight₁_le d B

/-! ## ★★組へ持ち上げる -/

/-- ★★**次数有界・高さ有界の代数的数の組は有限個**。

★`Set.Finite.pi` そのものである——1 成分の有限性(`boundedAlg_finite`)から直ちに出る。 -/
theorem finite_pi_boundedAlg {ι : Type*} [Finite ι] (d : ℕ) (B : ℝ) :
    {y : ι → ℂ | ∀ i, y i ∈ boundedAlg d B}.Finite := by
  have heq : {y : ι → ℂ | ∀ i, y i ∈ boundedAlg d B}
      = Set.pi Set.univ (fun _ : ι => boundedAlg d B) := by
    ext y; simp
  rw [heq]
  exact Set.Finite.pi (fun _ => boundedAlg_finite d B)

/-! ## ★★★★射影モデル経由の受け皿 -/

/-- ★★★★**次数有界・高さ有界の組へ単射に写る集合は有限**。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

★これが `Proposition 1.4, (iv)` を射影モデル経由で出すときの受け皿である——
`f` が「点をその射影座標に送る写像」、`hinj` が「埋め込みが単射」、
`hf` が「高さが `B` で抑えられているなら各座標の高さも抑えられる」に対応する。 -/
theorem finite_of_injOn_boundedAlg {α : Type*} {ι : Type*} [Finite ι] {d : ℕ} {B : ℝ}
    (f : α → (ι → ℂ)) (S : Set α) (hinj : Set.InjOn f S)
    (hf : ∀ x ∈ S, ∀ i, f x i ∈ boundedAlg d B) : S.Finite := by
  have himg : (f '' S).Finite := by
    refine Set.Finite.subset (finite_pi_boundedAlg (ι := ι) d B) ?_
    rintro y ⟨x, hx, rfl⟩
    exact fun i => hf x hx i
  exact Set.Finite.of_finite_image himg hinj

/-! ## ★出典の紐付け(`.src`) -/

def boundedAlg_finite.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(古典的 Northcott——1 成分)",
    sectionId := "genell-prop-1-4" }

def finite_pi_boundedAlg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(ℙⁿ の場合に要る組版)",
    sectionId := "genell-prop-1-4" }

def finite_of_injOn_boundedAlg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(射影埋め込み経由で X 上の Northcott を出す受け皿)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
