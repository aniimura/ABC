import ABC3.Found.GaloisRep.TatePhiHom

/-!
# Galois (G6) 第 293 ブロック —— **★★★★★★★★★★準同型性の仮定が類の言葉になった**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★到達点

> `c₁c₂c₃ = 1`、正規化代表元の積が `Q`、3 類が互いに異なり単位元でない —— このとき
> **`Φ(c₁) + Φ(c₂) + Φ(c₃) = 0`**(`tatePhi_add_add_eq_zero'`)

★★★第 292 では `a(c₁)a(c₂)a(c₃) = q` を仮定していたが、それが
**正規化代表元の積が `Q`** という類の言葉になった。

## ★★★★★★★正規化代表元の積は `Q^n`(`n = 0, 1, 2`)

`c₁c₂c₃ = 1` なら `u₁u₂u₃ ∈ q^ℤ`、すなわち `u₁u₂u₃ = Q^n`。付値を取れば

    n·v(Q) = v(u₁) + v(u₂) + v(u₃) ∈ [0, 3·v(Q))

★★したがって **`n ∈ {0, 1, 2}`**(`normRep_prod_zpow`)。

| `n` | 意味 | 扱い |
|---|---|---|
| `1` | 母数の積がちょうど `q` | **本ブロック**(第 292 に渡す) |
| `2` | 逆元を取れば `n = 1` | `Φ(c⁻¹) = −Φ(c)`(第 291)で還元 |
| `0` | 3 つとも単元 | 単元 2 つの群法則(第 273) |

## ★★★★積の条件は「単位元でない」から出る

`c₁c₂ = 1` は `c₃ = 1` と同値である(`c₁c₂c₃ = 1` だから)。★したがって第 292 の
`c_i·c_j ≠ 1` という 3 つの仮定は、**`c_k ≠ 1` の 3 つに吸収される**。
★★仮定が「単位元でない」「互いに異なる」の 6 本だけになった。

## ★実例のダイヤモンド(`tools/lean-idioms.md` に追記)

`Kˣ ⧸ Subgroup.zpowers Q` では `rw [one_mul]` が当たらない(`MulOneClass` の実例が
2 経路で来る)。`exact (one_mul c).symm.trans h` と**項の水準**で書けば通る。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `normRep_prod_zpow` | ★★★★★★★正規化代表元の積は `Q^n`(`n = 0,1,2`) |
| `tateAOf_prod_eq` | ★★★★★★★★`n = 1` なら母数の積が `q` |
| `tatePhi_add_add_eq_zero'` | ★★★★★★★★★★**準同型性(`n = 1` の場合)** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real WeierstrassCurve WeierstrassCurve.Affine QuotientGroup

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [Algebra R K]

/-! ## ★★★★★★★正規化代表元の積は `Q^n` -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**正規化代表元の積は `Q^n`(`n = 0, 1, 2`)**。

★付値の和が `[0, 3v(Q))` に入ることから `n` の範囲が出る。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem normRep_prod_zpow (S : TateSetup R I K) (c₁ c₂ c₃ : Kˣ ⧸ Subgroup.zpowers S.Q)
    (hprod : c₁ * c₂ * c₃ = 1) :
    ∃ n : ℤ, 0 ≤ n ∧ n ≤ 2 ∧
      normRep S.v S.Q S.hQ c₁ * normRep S.v S.Q S.hQ c₂ * normRep S.v S.Q S.hQ c₃ = S.Q ^ n := by
  set u₁ := normRep S.v S.Q S.hQ c₁ with hu1
  set u₂ := normRep S.v S.Q S.hQ c₂ with hu2
  set u₃ := normRep S.v S.Q S.hQ c₃ with hu3
  have hmk : (QuotientGroup.mk (u₁ * u₂ * u₃) : Kˣ ⧸ Subgroup.zpowers S.Q) = 1 := by
    rw [QuotientGroup.mk_mul, QuotientGroup.mk_mul, normRep_mk, normRep_mk, normRep_mk]
    exact hprod
  obtain ⟨n, hn⟩ := (QuotientGroup.eq_one_iff _).1 hmk
  have hn' : S.Q ^ n = u₁ * u₂ * u₃ := hn
  have hv : vAdd S.v u₁ + vAdd S.v u₂ + vAdd S.v u₃ = n * vAdd S.v S.Q := by
    rw [← vAdd_mul, ← vAdd_mul, ← hn', vAdd_zpow]
  have h1 := normRep_nonneg S.v S.Q S.hQ c₁
  have h2 := normRep_nonneg S.v S.Q S.hQ c₂
  have h3 := normRep_nonneg S.v S.Q S.hQ c₃
  have g1 := normRep_lt S.v S.Q S.hQ c₁
  have g2 := normRep_lt S.v S.Q S.hQ c₂
  have g3 := normRep_lt S.v S.Q S.hQ c₃
  have hQ := S.hQ
  refine ⟨n, ?_, ?_, hn'.symm⟩
  · nlinarith [hv, h1, h2, h3, hQ]
  · nlinarith [hv, g1, g2, g3, hQ]

/-- ★★★★★★★★**`n = 1` の場合は母数の積が `q` になる**。 -/
theorem tateAOf_prod_eq (S : TateSetup R I K) (c₁ c₂ c₃ : Kˣ ⧸ Subgroup.zpowers S.Q)
    (h1 : normRep S.v S.Q S.hQ c₁ * normRep S.v S.Q S.hQ c₂ * normRep S.v S.Q S.hQ c₃ = S.Q) :
    tateAOf S c₁ * tateAOf S c₂ * tateAOf S c₃ = S.q := by
  refine S.hinj ?_
  rw [map_mul, map_mul, (tateAOf_spec S c₁).1, (tateAOf_spec S c₂).1, (tateAOf_spec S c₃).1,
    S.hQq, ← Units.val_mul, ← Units.val_mul, h1]

/-! ## ★★★★★★★★★★準同型性(`n = 1` の場合) -/

section Hom

variable [DecidableEq K]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**準同型性(`n = 1` の場合)**——仮定は「単位元でない」「互いに異なる」。

★`c_i·c_j = 1` は `c_k = 1` と同値なので、第 292 の 3 つの積の仮定は吸収される。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tatePhi_add_add_eq_zero' (S : TateSetup R I K) (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (c₁ c₂ c₃ : Kˣ ⧸ Subgroup.zpowers S.Q) (hprod : c₁ * c₂ * c₃ = 1)
    (h1 : normRep S.v S.Q S.hQ c₁ * normRep S.v S.Q S.hQ c₂ * normRep S.v S.Q S.hQ c₃ = S.Q)
    (hn1 : c₁ ≠ 1) (hn2 : c₂ ≠ 1) (hn3 : c₃ ≠ 1)
    (h12 : c₁ ≠ c₂) (h13 : c₁ ≠ c₃) (h23 : c₂ ≠ c₃) :
    tatePhi S hΔ c₁ + tatePhi S hΔ c₂ + tatePhi S hΔ c₃ = 0 := by
  refine tatePhi_add_add_eq_zero S hloc hvR hvI hq0 hΔ c₁ c₂ c₃
    (tateAOf_prod_eq S c₁ c₂ c₃ h1) hn1 hn2 hn3 h12 h13 h23 ?_ ?_ ?_
  · intro h
    refine hn3 ?_
    have h2 := hprod
    rw [h] at h2
    exact (one_mul c₃).symm.trans h2
  · intro h
    refine hn2 ?_
    have h2 : c₂ * (c₁ * c₃) = 1 := by
      rw [← hprod, mul_comm c₁ c₂, mul_assoc]
    rw [h] at h2
    exact (mul_one c₂).symm.trans h2
  · intro h
    refine hn1 ?_
    have h2 : c₁ * (c₂ * c₃) = 1 := by
      rw [← hprod, mul_assoc]
    rw [h] at h2
    exact (mul_one c₁).symm.trans h2

end Hom

/-! ## ★出典の紐付け(`.src`) -/

def normRep_prod_zpow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——正規化代表元の積は Q^n)",
    sectionId := "genell-def-3-3" }

def tatePhi_add_add_eq_zero'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——準同型性(n = 1 の場合))",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
