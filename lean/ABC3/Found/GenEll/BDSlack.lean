import ABC3.Found.GenEll.BDClass

/-!
# [GenEll] 「`≤` を `≲` に変えるもの」(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.10。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

## ★★★★★★★★★これが `Proposition 1.7` の証明の骨である

原文 p.10 の証明は**明示的に 2 つに分けている**(2026-08-16 の 260 dpi 目視で確認):

- prime-to-`Σ` の部分は **`=` と `≤` であって `≲` ではない**
- `Σ` の上の部分は **`≥ 0` であって `≳` ではない**

★すなわち**厳密な不等式が `Σ` の外で成り立ち、`Σ` の上の食い違いが一様に有界**である。
★★**その 2 つから `≲`(= `BDle`)が出る**——それが本ファイルの補題である。

★★★`Σ` 上の食い違いが `Σ_{q∈Σ} log q` で一様に抑えられることは
`SigmaBound.lean` / `LogCondSigma.lean` で取ってある(第 412–415)。
**点にも定義体にも依らない定数**であることが要点である。
-/

namespace ABC3.Found.GenEll

variable {P : Type*}

/-- ★★★★★★★★★**厳密な不等式 + 一様に有界な食い違い ⟹ `BDle`**。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

`β x ≤ α x + slack x` が各点で成り立ち、`slack` が定数 `C` で抑えられるなら `α ≲ β`。

★**`BDle` の定数はそのまま `C` である**——`Σ` の上の寄与の限界がそのまま
`≲` の定数になる。これが「有限個の素数を除けば」という緩みが
BD-class に吸収される、ということの中身である。 -/
theorem bdle_of_bounded_slack (α β slack : P → ℝ) (C : ℝ)
    (h : ∀ x, β x ≤ α x + slack x) (hC : ∀ x, slack x ≤ C) :
    BDle α β :=
  ⟨C, fun x => by have h1 := h x; have h2 := hC x; linarith⟩

/-- ★**同じことを `BDge` の側で**。 -/
theorem bdge_of_bounded_slack (α β slack : P → ℝ) (C : ℝ)
    (h : ∀ x, α x ≤ β x + slack x) (hC : ∀ x, slack x ≤ C) :
    BDge α β :=
  ⟨C, fun x => by have h1 := h x; have h2 := hC x; linarith⟩

/-! ## ★出典の紐付け(`.src`) -/

def bdle_of_bounded_slack.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i)(prime-to-Σ の厳密な不等式と Σ 上の有界性から ≲ が出る)",
    sectionId := "genell-prop-1-7" }

end ABC3.Found.GenEll
