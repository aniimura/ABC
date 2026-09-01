/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Thm38ZetaPi
import ABC3.Meta.Claim

/-!
# 第 1277 ブロック —— **`π` は捩れから作れる（体を拡げなくてよい）**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★★★★★これは何か——第 1274 の道筋の**訂正**

第 1274 では「`K ≔ L_v(ζ_l, q^{1/l})` を建てる」（完備 DVR の完全分岐拡大）を
残る節点 B1 に挙げた。★しかし**それは要らない**。

☆`Lemma 3.2, (i)` の証明（在庫）が使っている Bezout の技がそのまま効く:

> `x^l = Q^k` かつ `l ∤ k` なら、`π ≔ x^b · Q^a`（`a l + b k = 1`）が `π^l = Q` を満たす。

★★★すなわち **`π` は「`Q` の `l` 乗根を添加した体」ではなく、
`l`-捩れ点の座標が住む体の中に既にある**。

☆そして `l`-捩れの類が `μ_l` の像に収まらなければ（個数 `l²` > `l` から出る）、
`l ∤ k` なる `x` が必ずある。

★★★したがって道筋は
「**大域の数体を `L′ ≔ L(ζ_l, E[l])` に取り替え、その完備化で `mkTateSetup` を使う**」
となり、局所体の建設（Eisenstein 拡大）は**まったく要らない**。
-/

namespace ABC3.Found.GenEll

open ABC3.Meta

/-- ★★★★★★★★★★★★
**`x^l = Q^k` かつ `l ∤ k` なら `Q` は `l` 乗になる**——★**無条件**（第 1277）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`lemma_3_2_i`（在庫）の証明の中にある Bezout の技を取り出したものである。 -/
theorem exists_pow_eq_of_not_dvd {G : Type*} [CommGroup G] {l : ℕ} (hl : l.Prime)
    (Q x : G) (k : ℤ) (hx : x ^ l = Q ^ k) (hk : ¬ ((l : ℤ) ∣ k)) :
    ∃ π : G, π ^ l = Q := by
  have hp : Prime (l : ℤ) := Nat.prime_iff_prime_int.mp hl
  obtain ⟨a, b, hab⟩ := (Prime.coprime_iff_not_dvd hp).mpr hk
  refine ⟨x ^ b * Q ^ a, ?_⟩
  have hxl : (x ^ b) ^ (l : ℕ) = Q ^ (k * b) := by
    rw [← zpow_natCast (x ^ b) l, ← zpow_mul, mul_comm b (l : ℤ), zpow_mul, zpow_natCast, hx,
      ← zpow_mul]
  rw [mul_pow, hxl, ← zpow_natCast (Q ^ a) l, ← zpow_mul, ← zpow_add]
  have he : k * b + a * (l : ℤ) = 1 := by linarith [hab]
  rw [he, zpow_one]

/-- ★★★★★★★★★★★★★★★★
**`μ_l` に収まらない `l`-捩れの類があれば `Q` は `l` 乗になる**——★**無条件**（第 1277）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`c^l = 1` なら `x^l = Q^k`。★`l ∣ k` なら `x·Q^{-k/l}` は `l` 乗根なので `ζ^a`、
すなわち `c = [ζ^a]` になる。☆そうでないと仮定したので `l ∤ k` である。

★★★これが「体を拡げずに `π` を作る」段である。 -/
theorem exists_pow_eq_of_torsion_not_mu {G : Type*} [CommGroup G] {l : ℕ} (hl : l.Prime)
    (Q ζ : G) (hμ : ∀ y : G, y ^ l = 1 → ∃ a : ℤ, y = ζ ^ a)
    (c : G ⧸ Subgroup.zpowers Q) (hc : c ^ l = 1)
    (hcnot : ∀ a : ℤ, c ≠ QuotientGroup.mk (ζ ^ a)) :
    ∃ π : G, π ^ l = Q := by
  obtain ⟨x, rfl⟩ := QuotientGroup.mk_surjective c
  obtain ⟨k, hk⟩ := (quotient_pow_eq_one_iff Q l x).mp hc
  by_cases hdvd : (l : ℤ) ∣ k
  · exfalso
    obtain ⟨m, rfl⟩ := hdvd
    have hroot : (x * Q ^ (-m)) ^ l = 1 := by
      rw [mul_pow, hk, ← zpow_natCast (Q ^ (-m)) l, ← zpow_mul, ← zpow_add]
      have he : (l : ℤ) * m + -m * (l : ℤ) = 0 := by ring
      rw [he, zpow_zero]
    obtain ⟨a, ha⟩ := hμ _ hroot
    refine hcnot a ?_
    have hx : x = ζ ^ a * Q ^ m := by
      have h2 := congrArg (fun z : G => z * Q ^ m) ha
      rw [mul_assoc, ← zpow_add] at h2
      simpa using h2
    rw [hx]
    refine QuotientGroup.eq.2 ⟨-m, ?_⟩
    group
  · exact exists_pow_eq_of_not_dvd hl Q x k hk hdvd

/-! ## ★出典の紐付け(`.src`) -/

def exists_pow_eq_of_not_dvd.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(x^l = Q^k かつ l ∤ k なら Q は l 乗。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_pow_eq_of_torsion_not_mu.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(μ_l に収まらない l-捩れの類があれば Q は l 乗。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_pow_eq_of_torsion_not_mu.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1277）**——第 1274 の道筋の**訂正**である。" ++
       "☆`π` は「`Q` の `l` 乗根を添加した体」ではなく、" ++
       "**`l`-捩れ点の座標が住む体の中に既にある**（Bezout）。" ++
       "★★★したがって局所体の建設（完備 DVR の Eisenstein 拡大、第 1274 の B1）は" ++
       "**まったく要らない**——大域の数体を `L′ ≔ L(ζ_l, E[l])` に取り替え、" ++
       "その完備化で `mkTateSetup`（在庫）を使えばよい。") 3 ]

end ABC3.Found.GenEll
