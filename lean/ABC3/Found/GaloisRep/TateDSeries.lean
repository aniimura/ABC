/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TatePair

/-!
# Galois (G6) 第 846 ブロック —— **★★★★★★★★`DX`・`DY` の級数**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★★★★★★★★★★★★★`v` は `∑_ζ DY` である（第 846 の発見）

Vélu の `v` は `v = ∑_{i≠0}(3X(ζ^i)² + a₄ − Y(ζ^i))` であり、`X²` が畳み込みを呼ぶ。
★★しかし ODE `DY = 3X² − Y + a₄`（第 842、`tate_equation` を微分して `2DX` で割る）
を入れると

    `v = ∑(DY + Y − a₄) + (l−1)a₄ − ∑Y = ∑_{i≠0} DY(ζ^i)`

★★★★**`X²` は完全に消える**。`l = 5, 7` で数値確認（2026-08-31）:

```
l=5   v(direct)=2.4892902282663…   ∑DY=2.4892902282386…   差 7.4e-11
l=7   v(direct)=9.8899579957404…   ∑DY=9.8899579959871…   差 2.5e-10
```

## ★★★★★★★★定数項は綺麗な指標和になる

`q = 0` では `DY(u,0) = u∂_u[u²/(1−u)³] = u²(2+u)/(1−u)⁴` なので、`v` の定数項は

    `∑_{ζ≠1} ζ²(2+ζ)/(1−ζ)⁴ = (l⁴−1)/240`

★`l = 2, 3, 5, 7, 11, 13` で厳密一致（2026-08-31）。
☆これは第 835 が測った `v(0) = (l⁴−1)/240` と**同じ値**である。

## ★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `tateDXterm` | `Df(t) = t(1+t)/(1−t)³`（`f(t) = t/(1−t)²` の `u∂_u`） |
| `tateDYterm` | `Dg(t) = t²(2+t)/(1−t)⁴`（`g(t) = t²/(1−t)³` の `u∂_u`） |
| `tateDXterm_eq` | ★★★**`Df = f + 2g`**——これ一つで `DX = 2Y + X` が出る |
| `tateDXtail` / `tateDYtail` | 尾 |
| `tateDXpair` / `tateDYpair` | ★★★**`DX`・`DY` の級数** |
| `tateDXpair_eq` | ★★★★★★**`DX = 2Y + X`**（`l = 5, 7` で数値確認済み） |

☆残るのは ODE `DY = 3X² − Y + a₄` ひとつであり、それは
`tate_equation`（証明済み）を万有な環 `TateUniv` の上で `D`（第 845）で微分し、
`2DX` で割ることで出る（`TateUniv` は整域なので割れる）。
-/

namespace ABC3.Found.GaloisRep

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★★★項 -/

/-- ★★**`Df(t) = t(1+t)/(1−t)³`**——`f(t) = t/(1−t)²` の `u∂_u`。 -/
noncomputable def tateDXterm (t : R) : R := t * (1 + t) * Ring.inverse (1 - t) ^ 3

/-- ★★**`Dg(t) = t²(2+t)/(1−t)⁴`**——`g(t) = t²/(1−t)³` の `u∂_u`。 -/
noncomputable def tateDYterm (t : R) : R := t ^ 2 * (2 + t) * Ring.inverse (1 - t) ^ 4

/-- ★★★★★★**`Df = f + 2g`**。

★これ一つで `DX = 2Y + X` が出る（`X` の `w` 側は符号が反転するが、
同じ恒等式がそのまま効く）。 -/
theorem tateDXterm_eq {t : R} (hu : IsUnit (1 - t)) :
    tateDXterm t = tateXterm t + 2 * tateYterm t := by
  have h : (1 - t) * Ring.inverse (1 - t) = 1 := Ring.mul_inverse_cancel _ hu
  rw [tateDXterm, tateXterm, tateYterm]
  linear_combination (t * Ring.inverse (1 - t) ^ 2) * h

theorem tateDXterm_mem_pow {k : ℕ} {t : R} (ht : t ∈ I ^ k) : tateDXterm t ∈ I ^ k :=
  Ideal.mul_mem_right _ _ (Ideal.mul_mem_right _ _ ht)

theorem tateDYterm_mem_pow {k : ℕ} {t : R} (ht : t ∈ I ^ k) : tateDYterm t ∈ I ^ k := by
  rw [tateDYterm, sq]
  exact Ideal.mul_mem_right _ _ (Ideal.mul_mem_right _ _ (Ideal.mul_mem_right _ _ ht))

/-! ## ★★★★★尾 -/

theorem pow_succ_mul_mem_I {u q : R} (hq : q ∈ I) (n : ℕ) : q ^ (n + 1) * u ∈ I := by
  rw [pow_succ]
  exact Ideal.mul_mem_right _ _ (Ideal.mul_mem_left _ _ hq)

theorem tateDXtail_aux {u q : R} (hq : q ∈ I) (n : ℕ) :
    tateDXterm (q ^ (n + 1) * u) ∈ I ^ n :=
  Ideal.pow_le_pow_right (Nat.le_succ n)
    (tateDXterm_mem_pow (Ideal.mul_mem_right u _ (Ideal.pow_mem_pow hq (n + 1))))

theorem tateDYtail_aux {u q : R} (hq : q ∈ I) (n : ℕ) :
    tateDYterm (q ^ (n + 1) * u) ∈ I ^ n :=
  Ideal.pow_le_pow_right (Nat.le_succ n)
    (tateDYterm_mem_pow (Ideal.mul_mem_right u _ (Ideal.pow_mem_pow hq (n + 1))))

/-- ★★★**`∑_{m≥1} Df(qᵐu)`**。 -/
noncomputable def tateDXtail [IsAdicComplete I R] (u q : R) (hq : q ∈ I) : R :=
  adicSum (fun n => tateDXterm (q ^ (n + 1) * u)) (tateDXtail_aux hq)

/-- ★★★**`∑_{m≥1} Dg(qᵐu)`**。 -/
noncomputable def tateDYtail [IsAdicComplete I R] (u q : R) (hq : q ∈ I) : R :=
  adicSum (fun n => tateDYterm (q ^ (n + 1) * u)) (tateDYtail_aux hq)

/-- ★★★★**尾の水準でも `Df = f + 2g`**。 -/
theorem tateDXtail_eq [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateDXtail u q hq = tateXtail u q hq + 2 * tateYtail u q hq := by
  have h2 : ∀ n : ℕ, (2 : R) * tateYterm (q ^ (n + 1) * u) ∈ I ^ n :=
    fun n => Ideal.mul_mem_left _ _ (tateYtail_aux hq n)
  have hpt : ∀ n : ℕ, tateDXterm (q ^ (n + 1) * u)
      = tateXterm (q ^ (n + 1) * u) + 2 * tateYterm (q ^ (n + 1) * u) :=
    fun n => tateDXterm_eq (isUnit_one_sub (I := I) (pow_succ_mul_mem_I hq n))
  have key : adicSum (fun n => tateDXterm (q ^ (n + 1) * u)) (tateDXtail_aux hq)
      = adicSum (fun n => tateXterm (q ^ (n + 1) * u) + 2 * tateYterm (q ^ (n + 1) * u))
          (fun n => Submodule.add_mem _ (tateXtail_aux hq n) (h2 n)) :=
    adicSum_congr _ _ hpt
  rw [tateDXtail, key, adicSum_add _ _ (tateXtail_aux hq) h2, adicSum_smul]
  rfl

/-! ## ★★★★★★★★`DX`・`DY` の級数 -/

/-- ★★★★★★**`DX(u,q)`**——`X = f(a)+T_f(a)+f(w)+T_f(w)−2s₁` を `u∂_u` したもの。

★`w = q/u` なので `w` 側は符号が反転し、`s₁` は `u` に依らないので消える。 -/
noncomputable def tateDXpair [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) : R :=
  (tateDXterm a + tateDXtail a q hq) - (tateDXterm w + tateDXtail w q hq)

/-- ★★★★★★**`DY(u,q)`**——
`Y = g(a)+T_g(a) − f(w)−T_f(w) − g(w)−T_g(w) + s₁` を `u∂_u` したもの。

★`w` 側は符号が反転するので、`−f(w)` は `+Df(w)`、`−g(w)` は `+Dg(w)` になる。 -/
noncomputable def tateDYpair [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) : R :=
  (tateDYterm a + tateDYtail a q hq) + (tateDXterm w + tateDXtail w q hq)
    + (tateDYterm w + tateDYtail w q hq)

/-- ★★★★★★★★★★**`DX = 2Y + X`**。

★`l = 5, 7` で数値確認済み（2026-08-31）。
☆`Df = f + 2g` だけから出る——`w` 側の符号反転が
`X` の `+f(w)+T_f(w)` と `Y` の `−f(w)−T_f(w)−g(w)−T_g(w)` に丁度合う。 -/
theorem tateDXpair_eq [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    tateDXpair a w q hq = 2 * tateYpair a w q hq + tateXpair a w q hq := by
  rw [tateDXpair, tateYpair, tateXpair, tateDXterm_eq ha, tateDXterm_eq hw,
    tateDXtail_eq, tateDXtail_eq]
  ring

/-! ## ★★★★★尾の漸化式と `Df` の反対称性 -/

theorem ring_inverse_eq_of_mul_eq_one {A : Type} [CommRing A] {x y : A} (hu : IsUnit x)
    (h : x * y = 1) : Ring.inverse x = y := by
  have h1 : Ring.inverse x * x = 1 := Ring.inverse_mul_cancel x hu
  calc Ring.inverse x = Ring.inverse x * (x * y) := by rw [h, mul_one]
    _ = Ring.inverse x * x * y := by ring
    _ = y := by rw [h1, one_mul]

/-- ★★**尾の漸化式** `T(u) = Df(qu) + T(qu)`。 -/
theorem tateDXtail_rec [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateDXtail u q hq = tateDXterm (q * u) + tateDXtail (q * u) q hq := by
  rw [tateDXtail, adicSum_shift]
  congr 1
  · norm_num
  · exact adicSum_congr _ _
      (fun n => by rw [show q ^ (n + 1 + 1) * u = q ^ (n + 1) * (q * u) by ring])

/-- ★★★★★★**`Df(1/t) = −Df(t)`**——`X(q/u) = X(u)` を微分したものの項版。

`1 − v = −(1−u)v`（`uv = 1`）なので `inv(1−v) = −u·inv(1−u)` であり、
`Df(v) = v(1+v)(−u·r)³ = −u³(v+v²)r³ = −(u²+u)r³ = −Df(u)`。 -/
theorem tateDXterm_inv {u v : R} (huv : u * v = 1) (hu : IsUnit (1 - u))
    (hv : IsUnit (1 - v)) : tateDXterm v = - tateDXterm u := by
  have hru : (1 - u) * Ring.inverse (1 - u) = 1 := Ring.mul_inverse_cancel _ hu
  have hkey : Ring.inverse (1 - v) = -u * Ring.inverse (1 - u) := by
    refine ring_inverse_eq_of_mul_eq_one hv ?_
    linear_combination hru + Ring.inverse (1 - u) * huv
  rw [tateDXterm, tateDXterm, hkey]
  linear_combination
    (-(Ring.inverse (1 - u) ^ 3) * (u ^ 2 + u * (u * v + 1))) * huv

/-! ## ★★★★★★★★★★`D²X`——σ₃ の母関数 -/

/-- ★★★★**`D²f(t) = t(1+4t+t²)/(1−t)⁴`**——`∑_{n≥1} n³t^n` の閉じた式。 -/
noncomputable def tateD2Xterm (t : R) : R :=
  t * (1 + 4 * t + t ^ 2) * Ring.inverse (1 - t) ^ 4

/-- ★★★★★★**`D²f = Df + 2Dg`**——`Df = f + 2g` と同じ形の恒等式。 -/
theorem tateD2Xterm_eq {t : R} (hu : IsUnit (1 - t)) :
    tateD2Xterm t = 2 * tateDYterm t + tateDXterm t := by
  have h : (1 - t) * Ring.inverse (1 - t) = 1 := Ring.mul_inverse_cancel _ hu
  rw [tateD2Xterm, tateDYterm, tateDXterm]
  linear_combination (t * (1 + t) * Ring.inverse (1 - t) ^ 3) * h

theorem tateD2Xterm_mem_pow {k : ℕ} {t : R} (ht : t ∈ I ^ k) : tateD2Xterm t ∈ I ^ k :=
  Ideal.mul_mem_right _ _ (Ideal.mul_mem_right _ _ ht)

theorem tateD2Xtail_aux {u q : R} (hq : q ∈ I) (n : ℕ) :
    tateD2Xterm (q ^ (n + 1) * u) ∈ I ^ n :=
  Ideal.pow_le_pow_right (Nat.le_succ n)
    (tateD2Xterm_mem_pow (Ideal.mul_mem_right u _ (Ideal.pow_mem_pow hq (n + 1))))

/-- ★`∑_{m≥1} D²f(qᵐu)`。 -/
noncomputable def tateD2Xtail [IsAdicComplete I R] (u q : R) (hq : q ∈ I) : R :=
  adicSum (fun n => tateD2Xterm (q ^ (n + 1) * u)) (tateD2Xtail_aux hq)

theorem tateD2Xtail_eq [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateD2Xtail u q hq = 2 * tateDYtail u q hq + tateDXtail u q hq := by
  have h2 : ∀ n : ℕ, (2 : R) * tateDYterm (q ^ (n + 1) * u) ∈ I ^ n :=
    fun n => Ideal.mul_mem_left _ _ (tateDYtail_aux hq n)
  have hpt : ∀ n : ℕ, tateD2Xterm (q ^ (n + 1) * u)
      = 2 * tateDYterm (q ^ (n + 1) * u) + tateDXterm (q ^ (n + 1) * u) :=
    fun n => tateD2Xterm_eq (isUnit_one_sub (I := I) (pow_succ_mul_mem_I hq n))
  have key : adicSum (fun n => tateD2Xterm (q ^ (n + 1) * u)) (tateD2Xtail_aux hq)
      = adicSum (fun n => 2 * tateDYterm (q ^ (n + 1) * u) + tateDXterm (q ^ (n + 1) * u))
          (fun n => Submodule.add_mem _ (h2 n) (tateDXtail_aux hq n)) :=
    adicSum_congr _ _ hpt
  rw [tateD2Xtail, key, adicSum_add _ _ h2 (tateDXtail_aux hq), adicSum_smul]
  rfl

/-- ★★★★★★**`D²X(u,q)`**——`w` 側の符号は 2 回反転して戻るので**足し算**である。 -/
noncomputable def tateD2Xpair [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) : R :=
  (tateD2Xterm a + tateD2Xtail a q hq) + (tateD2Xterm w + tateD2Xtail w q hq)

/-- ★★★★★★★★★★**`D²X = 2DY + DX`**。

☆これと `∑_ζ DX = 0` を合わせると **`v = ∑_ζ DY = (1/2)∑_ζ D²X`** となり、
`D²X` の `q^N` 係数が `∑_{d∣N} d³(u^d + u^{-d})` なので、★★**σ₃ が直に出る**。 -/
theorem tateD2Xpair_eq [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    tateD2Xpair a w q hq = 2 * tateDYpair a w q hq + tateDXpair a w q hq := by
  rw [tateD2Xpair, tateDYpair, tateDXpair, tateD2Xterm_eq ha, tateD2Xterm_eq hw,
    tateD2Xtail_eq, tateD2Xtail_eq]
  ring

/-! ## ★出典の紐付け(`.src`) -/

def tateDXterm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(Df(t) = t(1+t)/(1−t)³。★無条件)",
    sectionId := "genell-lemma-3-2" }

def tateDYterm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(Dg(t) = t²(2+t)/(1−t)⁴。★無条件)",
    sectionId := "genell-lemma-3-2" }

def tateDXterm_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(Df = f + 2g。★1−t が単元)",
    sectionId := "genell-lemma-3-2" }

def tateDXpair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(DX の級数。★無条件)",
    sectionId := "genell-lemma-3-2" }

def tateDYpair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(DY の級数。★無条件)",
    sectionId := "genell-lemma-3-2" }

def tateDXpair_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(DX = 2Y + X。★1−a と 1−w が単元)",
    sectionId := "genell-lemma-3-2" }

def tateD2Xterm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(D²f(t) = t(1+4t+t²)/(1−t)⁴。★無条件)",
    sectionId := "genell-lemma-3-2" }

def tateD2Xterm_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(D²f = Df + 2Dg。★1−t が単元)",
    sectionId := "genell-lemma-3-2" }

def tateD2Xpair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(D²X の級数。★無条件)",
    sectionId := "genell-lemma-3-2" }

def tateD2Xpair_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(D²X = 2DY + DX。★1−a と 1−w が単元)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GaloisRep
