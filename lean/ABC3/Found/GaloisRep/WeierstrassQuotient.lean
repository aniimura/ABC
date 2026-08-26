import ABC3.Found.GaloisRep.WeierstrassPoleForm

/-!
# Galois (G6) 第 245 ブロック —— **★★★★★★★格子点での正則延長**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★因数分解を「商の形」に持ち上げる

第 244 の `pole_cancel_factor` は多項式の恒等式だった。本ブロックはそれを
**除算のある形**に翻訳し、さらに `℘`・`℘'` に代入して

    (1/4)·((℘'(z) − c')/(℘(z) − c))² − ℘(z)

が**格子点 `l₀` のまわりで解析的な関数に一致する**ことを示す
(`exists_analytic_quotient_sub`)。ここで `c := ℘(w)`、`c' := ℘'(w)` は定数でよい。

`t := z − l₀`、`α := D z − c'`、`β := E z − c` として第 244 の形に合わせると

    (℘'(z) − c')/(℘(z) − c) = (−2 + t³α)/(t·(1 + t²β))        `quotient_pole_form`
    (1/4)·(…)² − 1/t²      = R/(4·(1 + t²β)²)                  `sq_quotient_sub_pole`

であり、`℘(z) = 1/t² + E z` を引くと `1/t²` が消えて

    g z := R z /(4·(1 + t²β)²) − E z

が残る。`R`・`E` は `l₀` で解析的、分母は `l₀` で `4 ≠ 0` なので `g` は解析的。

## ★★★これが 3 種類の悪い点をまとめて処理する

加法定理 `F(z) = ℘(z+w) + ℘(z) + ℘(w) − (1/4)((℘'(z)−℘'(w))/(℘(z)−℘(w)))²` の
見かけの極は `z ∈ L`、`z ∈ L − w`、`z ≡ w (mod L)` の 3 種類だが、
どれも**同じ商の形**に帰着する。本ブロックの `exists_analytic_quotient_sub` は
そのうち第 1 種(`z ∈ L`)をそのまま与え、残る 2 種も同じ補題の平行移動で出る。

## ★分母が消えないことは `∀ᶠ` で持つ

`1 + t²β` は `l₀` で `1` なので、**`l₀` の近傍で** `≠ 0` である。等式を全 `z` で
主張せず `∀ᶠ z in 𝓝 l₀` の形にすることで、余計な仮定を呼び出し側に押し付けずに済む。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `quotient_pole_form` | ★★★★商の形——分母から `t` が 1 つ出る |
| `sq_quotient_sub_pole` | ★★★★★**極の相殺(商の形)** |
| `exists_analytic_quotient_sub` | ★★★★★★★**格子点での正則延長** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real PeriodPair

/-! ## ★★★★商の形 -/

/-- ★★★★**商の形**——`(−2/t³ + α)/(1/t² + β)` の分母から `t` が 1 つ出る。 -/
theorem quotient_pole_form (t α β : ℂ) (ht : t ≠ 0) (hM : 1 + t ^ 2 * β ≠ 0) :
    (-2 / t ^ 3 + α) / (1 / t ^ 2 + β) = (-2 + t ^ 3 * α) / (t * (1 + t ^ 2 * β)) := by
  have ht2 : t ^ 2 ≠ 0 := pow_ne_zero 2 ht
  have ht3 : t ^ 3 ≠ 0 := pow_ne_zero 3 ht
  have hden : (1 / t ^ 2 + β) ≠ 0 := by
    intro h
    apply hM
    field_simp at h ⊢
    linear_combination h
  field_simp

/-- ★★★★★**極の相殺(商の形)**——`1/t²` を引くと分母から `t` が消える。

第 244 の `pole_cancel_factor` を除算のある形に翻訳したもの。 -/
theorem sq_quotient_sub_pole (t α β : ℂ) (ht : t ≠ 0) (hM : 1 + t ^ 2 * β ≠ 0) :
    (1 / 4) * ((-2 + t ^ 3 * α) / (t * (1 + t ^ 2 * β))) ^ 2 - 1 / t ^ 2
      = (-4 * t * α + t ^ 4 * α ^ 2 - 8 * β - 4 * t ^ 2 * β ^ 2) / (4 * (1 + t ^ 2 * β) ^ 2) := by
  have ht2 : t ^ 2 ≠ 0 := pow_ne_zero 2 ht
  field_simp
  ring

/-! ## ★★★★★★★格子点での正則延長 -/

set_option maxHeartbeats 800000 in
/-- ★★★★★★★**格子点での正則延長**——`(1/4)·((℘' − c')/(℘ − c))² − ℘` は
格子点 `l₀` のまわりで解析的な関数に一致する。

加法定理の証明で「見かけの極が相殺する」ことの中身。`c`・`c'` は定数でよいので、
`c := ℘(w)`、`c' := ℘'(w)` を入れればそのまま使える。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem exists_analytic_quotient_sub (L : PeriodPair) (l₀ : L.lattice) (c c' : ℂ) :
    ∃ g : ℂ → ℂ, AnalyticAt ℂ g (l₀ : ℂ) ∧
      ∀ᶠ z in nhds (l₀ : ℂ), z ≠ (l₀ : ℂ) →
        (1 / 4) * ((L.derivWeierstrassP z - c') / (L.weierstrassP z - c)) ^ 2
          - L.weierstrassP z = g z := by
  obtain ⟨E, D, hE, hD, hP, hP'⟩ := exists_pole_form L l₀
  set num : ℂ → ℂ := fun z =>
    -4 * (z - (l₀ : ℂ)) * (D z - c') + (z - (l₀ : ℂ)) ^ 4 * (D z - c') ^ 2
      - 8 * (E z - c) - 4 * (z - (l₀ : ℂ)) ^ 2 * (E z - c) ^ 2 with hnumdef
  set den : ℂ → ℂ := fun z => 4 * (1 + (z - (l₀ : ℂ)) ^ 2 * (E z - c)) ^ 2 with hdendef
  have hid : AnalyticAt ℂ (fun z : ℂ => z - (l₀ : ℂ)) (l₀ : ℂ) :=
    analyticAt_id.sub analyticAt_const
  have hnum : AnalyticAt ℂ num (l₀ : ℂ) := by
    simp only [hnumdef]
    exact ((((analyticAt_const.mul hid).mul (hD.sub analyticAt_const)).add
      ((hid.pow 4).mul ((hD.sub analyticAt_const).pow 2))).sub
      (analyticAt_const.mul (hE.sub analyticAt_const))).sub
      ((analyticAt_const.mul (hid.pow 2)).mul ((hE.sub analyticAt_const).pow 2))
  have hden : AnalyticAt ℂ den (l₀ : ℂ) := by
    simp only [hdendef]
    exact analyticAt_const.mul ((analyticAt_const.add
      ((hid.pow 2).mul (hE.sub analyticAt_const))).pow 2)
  have hden0 : den (l₀ : ℂ) ≠ 0 := by
    simp only [hdendef]
    norm_num
  refine ⟨fun z => num z / den z - E z, (hnum.div hden hden0).sub hE, ?_⟩
  have hMcont : ContinuousAt (fun z : ℂ => 1 + (z - (l₀ : ℂ)) ^ 2 * (E z - c)) (l₀ : ℂ) := by
    have := hE.continuousAt
    fun_prop
  have hMne : ∀ᶠ z in nhds (l₀ : ℂ), (1 + (z - (l₀ : ℂ)) ^ 2 * (E z - c)) ≠ 0 := by
    refine hMcont.eventually_ne ?_
    norm_num
  filter_upwards [hMne] with z hMz hz
  have ht : z - (l₀ : ℂ) ≠ 0 := sub_ne_zero.2 hz
  have h1 : L.derivWeierstrassP z - c' = -2 / (z - (l₀ : ℂ)) ^ 3 + (D z - c') := by
    rw [hP' z]
    ring
  have h2 : L.weierstrassP z - c = 1 / (z - (l₀ : ℂ)) ^ 2 + (E z - c) := by
    rw [hP z]
    ring
  rw [h1, h2, quotient_pole_form _ _ _ ht hMz]
  have h3 := sq_quotient_sub_pole (z - (l₀ : ℂ)) (D z - c') (E z - c) ht hMz
  have h4 : L.weierstrassP z = 1 / (z - (l₀ : ℂ)) ^ 2 + E z := hP z
  rw [h4]
  simp only [hnumdef, hdendef]
  linear_combination h3

/-! ## ★出典の紐付け(`.src`) -/

def sq_quotient_sub_pole.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——極の相殺を商の形で)",
    sectionId := "genell-def-3-3" }

def exists_analytic_quotient_sub.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——格子点での正則延長)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
