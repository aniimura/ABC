import ABC3.Found.GaloisRep.TateEquationMod
import ABC3.Found.GaloisRep.TateQExp

/-!
# Galois (G6) 第 221 ブロック —— **★★★★★★★Tate 級数は環準同型で移る**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★段 6 の道具

第 220 で解析側の Tate 方程式が出た。これを形式側(任意の `I` 進完備環 `R`)へ移すには
**「万有な環で 0 なら任意の環で 0」**が要る。そのための**特殊化の道具**が本ブロックである。

★連続な環準同型 `φ : R →+* R'`(`I^n` を `J^n` に送る)に対して、
Tate 級数のすべての部品が `φ` で移ることを確かめる。

| 部品 | 移る理由 |
|---|---|
| `Ring.inverse x` | ★単元の上でのみ——`x` が単元なら `φ x` も単元 |
| `adicSum a` | ★★`adicSum_unique`(極限の一意性)——`φ` が `I^n → J^n` だから部分和の合同が移る |
| `evalAdic f q` | ★★同じ(`partialEval` は `Int.cast` と `q^n` だけ) |
| `tateXtail` / `tateYtail` | ★`adicSum` + 各項 |
| `tateXpair` / `tateYpair` | ★`1 − a`・`1 − w` が単元であることが要る |
| `tateCurveAt` | ★★`evalAdicHom` の合成(`WeierstrassCurve.map_map`) |

## ★★★`Ring.inverse` が要注意である

★★`Ring.inverse` は環準同型と**一般には可換でない**——`x` が単元でないとき既定値 `0` を返し、
`φ x` が単元になっても `φ 0 = 0 ≠ Ring.inverse (φ x)` となる。
★**単元であるという仮定が要る**(`map_ring_inverse`)。これが `tateXpair` の
移送に `IsUnit (1 − a)`・`IsUnit (1 − w)` を要求する理由である。
第 212 ブロックで `Ideal.Quotient` を避けたのも同じ理由だった。

## ★★到達点

    φ (tateDefect a w q hq) = tateDefect (φ a) (φ w) (φ q) hq'

★`tateDefect` は Weierstrass 方程式の左辺 − 右辺である。これが 0 であることが**葉 (b)**。
★★本ブロックにより、**万有な環で `tateDefect = 0` を示せば任意の完備環で従う**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `map_ring_inverse` | ★単元の上で `Ring.inverse` は移る |
| `smodEq_iff_sub_mem` | ★`SModEq (I^n • ⊤)` は `I^n` の元であること |
| `map_adicSum` | ★★★★★**`I` 進和は移る** |
| `map_evalAdic` | ★★★★★**冪級数の値は移る** |
| `map_tateXpair` / `map_tateYpair` | ★★★★★★**Tate 級数は移る** |
| `map_tateCurveAt` | ★★★★★Tate 曲線は移る |
| `map_tateDefect` | ★★★★★★★**方程式の差は移る** |
-/

namespace ABC3.Found.GaloisRep

variable {R R' : Type} [CommRing R] [CommRing R'] {I : Ideal R} {J : Ideal R'}

/-! ## ★`Ring.inverse` と `SModEq` -/

/-- ★環準同型は `Ring.inverse` と可換である——**単元の上でのみ**。 -/
theorem map_ring_inverse (φ : R →+* R') {x : R} (hx : IsUnit x) :
    φ (Ring.inverse x) = Ring.inverse (φ x) := by
  have hux : IsUnit (φ x) := hx.map φ
  have h : φ (Ring.inverse x) * φ x = 1 := by
    rw [← map_mul, Ring.inverse_mul_cancel _ hx, map_one]
  calc φ (Ring.inverse x) = φ (Ring.inverse x) * (φ x * Ring.inverse (φ x)) := by
        rw [Ring.mul_inverse_cancel _ hux, mul_one]
    _ = (φ (Ring.inverse x) * φ x) * Ring.inverse (φ x) := by ring
    _ = Ring.inverse (φ x) := by rw [h, one_mul]

/-- ★`I^n • ⊤` を法とする合同は差が `I^n` の元であることに等しい。 -/
theorem smodEq_iff_sub_mem {n : ℕ} (x y : R) :
    x ≡ y [SMOD (I ^ n • (⊤ : Submodule R R))] ↔ x - y ∈ I ^ n := by
  rw [SModEq.sub_mem]
  constructor
  · intro h
    simpa [Ideal.smul_top_eq_map] using h
  · intro h
    simpa [Ideal.smul_top_eq_map] using h

/-! ## ★項と部分和 -/

/-- ★`f(t)` は環準同型で移る(`1 − t` が単元のとき)。 -/
theorem map_tateXterm (φ : R →+* R') {t : R} (ht : IsUnit (1 - t)) :
    φ (tateXterm t) = tateXterm (φ t) := by
  rw [tateXterm, tateXterm, map_mul, map_pow, map_ring_inverse φ ht, map_sub, map_one]

/-- ★`g(t)` は環準同型で移る。 -/
theorem map_tateYterm (φ : R →+* R') {t : R} (ht : IsUnit (1 - t)) :
    φ (tateYterm t) = tateYterm (φ t) := by
  rw [tateYterm, tateYterm, map_mul, map_pow, map_pow, map_ring_inverse φ ht, map_sub, map_one]

theorem map_partialSum (φ : R →+* R') (a : ℕ → R) (N : ℕ) :
    φ (partialSum a N) = partialSum (fun n => φ (a n)) N := by
  rw [partialSum, partialSum, map_sum]

theorem map_partialEval (φ : R →+* R') (f : PowerSeries ℤ) (q : R) (N : ℕ) :
    φ (partialEval f q N) = partialEval f (φ q) N := by
  rw [partialEval, partialEval, map_sum]
  refine Finset.sum_congr rfl fun n _ => ?_
  rw [map_mul, map_pow, map_intCast]

/-! ## ★★★★★`I` 進和と冪級数の値 -/

/-- ★★★★★**`adicSum` は連続な環準同型で移る**。

★極限の一意性(`adicSum_unique`)から出る——部分和は明らかに移り、
`φ` が `I^n` を `J^n` に送るので合同も移る。 -/
theorem map_adicSum [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
    (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n)
    (a : ℕ → R) (ha : ∀ n, a n ∈ I ^ n) :
    φ (adicSum a ha) = adicSum (fun n => φ (a n)) (fun n => hφ n _ (ha n)) := by
  refine (adicSum_unique _ _ _ ?_).symm
  intro n
  rw [smodEq_iff_sub_mem, ← map_partialSum φ a n, ← map_sub]
  exact hφ n _ ((smodEq_iff_sub_mem _ _).1 (adicSum_spec a ha n))

/-- ★★★★★**冪級数の値は連続な環準同型で移る**。 -/
theorem map_evalAdic [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
    (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n)
    (f : PowerSeries ℤ) (q : R) (hq : q ∈ I) :
    φ (evalAdic f q hq) = evalAdic f (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  refine (evalAdic_unique _ _ _ _ (fun n => ?_)).symm
  rw [smodEq_iff_sub_mem, ← map_partialEval φ f q n, ← map_sub]
  exact hφ n _ ((smodEq_iff_sub_mem _ _).1 (evalAdic_spec f q hq n))

/-! ## ★★★★★★Tate 級数 -/

/-- ★★★★**片側の尾は環準同型で移る**。 -/
theorem map_tateXtail [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
    (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n)
    (u q : R) (hq : q ∈ I) :
    φ (tateXtail u q hq)
      = tateXtail (φ u) (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  rw [tateXtail, map_adicSum φ hφ, tateXtail]
  refine adicSum_congr _ _ fun n => ?_
  rw [map_tateXterm φ (isUnit_one_sub (I := I) (pow_succ_mul_mem (u := u) hq n)),
    map_mul, map_pow]

theorem map_tateYtail [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
    (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n)
    (u q : R) (hq : q ∈ I) :
    φ (tateYtail u q hq)
      = tateYtail (φ u) (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  rw [tateYtail, map_adicSum φ hφ, tateYtail]
  refine adicSum_congr _ _ fun n => ?_
  rw [map_tateYterm φ (isUnit_one_sub (I := I) (pow_succ_mul_mem (u := u) hq n)),
    map_mul, map_pow]

/-- ★★★★★★**`X(u,q)` は環準同型で移る**。 -/
theorem map_tateXpair [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
    (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n)
    (a w q : R) (hq : q ∈ I) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    φ (tateXpair a w q hq)
      = tateXpair (φ a) (φ w) (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  rw [tateXpair, tateXpair, map_sub, map_add, map_add, map_add, map_mul,
    map_tateXterm φ ha, map_tateXterm φ hw,
    map_tateXtail φ hφ a q hq, map_tateXtail φ hφ w q hq,
    map_evalAdic φ hφ _ q hq, map_ofNat φ 2]

/-- ★★★★★★**`Y(u,q)` は環準同型で移る**。 -/
theorem map_tateYpair [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
    (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n)
    (a w q : R) (hq : q ∈ I) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    φ (tateYpair a w q hq)
      = tateYpair (φ a) (φ w) (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  rw [tateYpair, tateYpair, map_add, map_sub, map_sub, map_add, map_add, map_add,
    map_tateYterm φ ha, map_tateYterm φ hw, map_tateXterm φ hw,
    map_tateYtail φ hφ a q hq, map_tateYtail φ hφ w q hq, map_tateXtail φ hφ w q hq,
    map_evalAdic φ hφ _ q hq]

/-! ## ★★★★★Tate 曲線 -/

theorem evalAdicHom_comp [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
    (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n) (q : R) (hq : q ∈ I) :
    φ.comp (evalAdicHom q hq)
      = evalAdicHom (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  refine RingHom.ext fun f => ?_
  exact map_evalAdic φ hφ f q hq

/-- ★★★★★**Tate 曲線は環準同型で移る**。 -/
theorem map_tateCurveAt [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
    (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n) (q : R) (hq : q ∈ I) :
    (tateCurveAt q hq).map φ
      = tateCurveAt (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  rw [tateCurveAt, WeierstrassCurve.map_map, evalAdicHom_comp φ hφ q hq, tateCurveAt]

theorem map_tateCurveAt_a4 [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
    (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n) (q : R) (hq : q ∈ I) :
    φ ((tateCurveAt q hq).a₄)
      = (tateCurveAt (φ q) (by simpa using hφ 1 q (by simpa using hq))).a₄ := by
  rw [← map_tateCurveAt φ hφ q hq]
  rfl

theorem map_tateCurveAt_a6 [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
    (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n) (q : R) (hq : q ∈ I) :
    φ ((tateCurveAt q hq).a₆)
      = (tateCurveAt (φ q) (by simpa using hφ 1 q (by simpa using hq))).a₆ := by
  rw [← map_tateCurveAt φ hφ q hq]
  rfl

/-! ## ★★★★★★★方程式の差 -/

/-- ★**Weierstrass 方程式の差**——これが 0 であることが葉 (b)。 -/
noncomputable def tateDefect [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) : R :=
  tateYpair a w q hq ^ 2 + tateXpair a w q hq * tateYpair a w q hq
    - (tateXpair a w q hq ^ 3 + (tateCurveAt q hq).a₄ * tateXpair a w q hq
      + (tateCurveAt q hq).a₆)

/-- ★★第 212 ブロックの帰結——差は `I` の元である。 -/
theorem tateDefect_mem [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (hw : w ∈ I)
    (hu : IsUnit (1 - a)) : tateDefect a w q hq ∈ I :=
  tate_equation_mem a w q hq hw hu

/-- ★★★★★★★**方程式の差は環準同型で移る**——これが特殊化の道具である。

★★これで「万有な環で 0」から「任意の完備環で 0」が出る。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem map_tateDefect [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
    (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n)
    (a w q : R) (hq : q ∈ I) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    φ (tateDefect a w q hq)
      = tateDefect (φ a) (φ w) (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  rw [tateDefect, tateDefect, map_sub, map_add, map_add, map_add, map_mul, map_mul,
    map_pow, map_pow, map_tateXpair φ hφ a w q hq ha hw,
    map_tateYpair φ hφ a w q hq ha hw,
    map_tateCurveAt_a4 φ hφ q hq, map_tateCurveAt_a6 φ hφ q hq]

/-- ★★★★★★**差が 0 であることは環準同型で移る**。 -/
theorem tateDefect_eq_zero_map [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
    (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n)
    (a w q : R) (hq : q ∈ I) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w))
    (h : tateDefect a w q hq = 0) :
    tateDefect (φ a) (φ w) (φ q) (by simpa using hφ 1 q (by simpa using hq)) = 0 := by
  rw [← map_tateDefect φ hφ a w q hq ha hw, h, map_zero]

/-- ★★★差が 0 なら Weierstrass 方程式が成り立つ。 -/
theorem tate_equation_of_defect [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (h : tateDefect a w q hq = 0) :
    tateYpair a w q hq ^ 2 + tateXpair a w q hq * tateYpair a w q hq
      = tateXpair a w q hq ^ 3 + (tateCurveAt q hq).a₄ * tateXpair a w q hq
        + (tateCurveAt q hq).a₆ := by
  rw [tateDefect, sub_eq_zero] at h
  exact h

/-! ## ★出典の紐付け(`.src`) -/

def map_adicSum.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——I 進和の函手性)",
    sectionId := "genell-def-3-3" }

def map_tateXpair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——Tate 級数の函手性)",
    sectionId := "genell-def-3-3" }

def map_tateDefect.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——方程式の差の函手性、特殊化の道具)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
