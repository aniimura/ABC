import ABC3.Found.GaloisRep.CrtStep

/-!
# Galois (G1) 第 71 ブロック —— **★★★★★★★★有限アーベル群の構造定理**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★これで群論の段が閉じる

    A 有限アーベル、n•x = 0、#A[d] = d² (d ∣ n)  ⟹  A ≃+ (ℤ/n)²

★`n` についての強い帰納法:

| 段 | 内容 |
|---|---|
| `n = 1` | ★`A` は自明 |
| `n > 1` | ★★`p = n.minFac`、`k = v_p(n)`、`m = n/p^k` |
| | ★★★`A ≃+ A[p^k] × A[m]`(第 70、中国剰余) |
| | ★★★★`A[p^k] ≃+ (ℤ/p^k)²`(第 69、素数冪) |
| | ★★★★★`A[m] ≃+ (ℤ/m)²`(帰納法、`m < n`) |
| | ★★★★★★`(ℤ/p^k)² × (ℤ/m)² ≃+ (ℤ/n)²`(第 70) |

★★★部分群の中の捩れが元の捩れと同じ(`kerSubEquiv`)なので、
仮定 `#A[d] = d²` が `A[p^k]` と `A[m]` にそのまま降りる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `addEquiv_zmod_sq` | ★★★★★★★★**`A ≃+ (ℤ/n)²`** |
-/

namespace ABC3.Found.GaloisRep

universe u

/-- ★★★★★★★★**有限アーベル群の構造定理**——
`n•x = 0` かつ `#A[d] = d²` (`d ∣ n`) なら `A ≃+ (ℤ/n)²`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★これが `Theorem 3.8` の群論的な核である。 -/
theorem addEquiv_zmod_sq {A : Type u} [AddCommGroup A] [Finite A] (n : ℕ) (hn : 1 ≤ n)
    (hexp : ∀ x : A, n • x = 0)
    (hcard : ∀ d : ℕ, d ∣ n → Nat.card (nsmulHom A d).ker = d ^ 2) :
    Nonempty ((ZMod n × ZMod n) ≃+ A) := by
  induction n using Nat.strong_induction_on generalizing A with
  | _ n ih =>
    rcases eq_or_lt_of_le hn with h1 | h1
    · have hn1 : n = 1 := h1.symm
      subst hn1
      have hsub : Subsingleton A := ⟨fun x y => by
        have hx := hexp x; have hy := hexp y; simp only [one_smul] at hx hy; rw [hx, hy]⟩
      refine ⟨AddEquiv.ofBijective (0 : (ZMod 1 × ZMod 1) →+ A) ⟨?_, ?_⟩⟩
      · intro x y _; exact Subsingleton.elim x y
      · intro y; exact ⟨0, Subsingleton.elim _ _⟩
    · have hn0 : n ≠ 0 := by omega
      set p := n.minFac with hp
      have hpp : p.Prime := Nat.minFac_prime (by omega)
      set k := n.factorization p with hk
      set m := n / p ^ k with hm
      have hsplit : p ^ k * m = n := Nat.ordProj_mul_ordCompl_eq_self n p
      have hco : Nat.Coprime (p ^ k) m := Nat.Coprime.pow_left _ (Nat.coprime_ordCompl hpp hn0)
      have hk1 : 1 ≤ k := Nat.Prime.factorization_pos_of_dvd hpp hn0 (Nat.minFac_dvd n)
      have hpkdvd : p ^ k ∣ n := ⟨m, hsplit.symm⟩
      have hmdvd : m ∣ n := ⟨p ^ k, by rw [← hsplit]; ring⟩
      have hpk2 : 2 ≤ p ^ k := le_trans hpp.two_le (Nat.le_self_pow (by omega) p)
      have hmpos : 0 < m := Nat.pos_of_ne_zero (fun h0 => by rw [h0, mul_zero] at hsplit; omega)
      have hmlt : m < n := by nlinarith [hsplit]
      obtain ⟨eCrt⟩ := addEquiv_crt (A := A) (p ^ k) m hco (by rw [hsplit]; exact hexp)
      have hexpP : ∀ x : (nsmulHom A (p ^ k)).ker, (p ^ k) • x = 0 := fun x => by ext; exact x.2
      have hcardP : Nat.card (nsmulHom A (p ^ k)).ker = (p ^ k) ^ 2 := hcard _ hpkdvd
      have hcardPk1 : Nat.card (nsmulHom ((nsmulHom A (p ^ k)).ker) (p ^ (k - 1))).ker
          = (p ^ (k - 1)) ^ 2 := by
        rw [Nat.card_congr (kerSubEquiv (p ^ k) (p ^ (k - 1)) (pow_dvd_pow p (by omega))).toEquiv]
        exact hcard _ (dvd_trans (pow_dvd_pow p (by omega)) hpkdvd)
      have hcardP1 : Nat.card (nsmulHom ((nsmulHom A (p ^ k)).ker) p).ker = p ^ 2 := by
        rw [Nat.card_congr (kerSubEquiv (p ^ k) p (dvd_pow_self p (by omega))).toEquiv]
        exact hcard _ (dvd_trans (dvd_pow_self p (by omega)) hpkdvd)
      obtain ⟨eP⟩ := addEquiv_prime_pow p k hpp hk1 hexpP hcardP hcardPk1 hcardP1
      have hexpM : ∀ x : (nsmulHom A m).ker, m • x = 0 := fun x => by ext; exact x.2
      have hcardM : ∀ d : ℕ, d ∣ m → Nat.card (nsmulHom ((nsmulHom A m).ker) d).ker = d ^ 2 := by
        intro d hd
        rw [Nat.card_congr (kerSubEquiv m d hd).toEquiv]
        exact hcard d (dvd_trans hd hmdvd)
      obtain ⟨eM⟩ := ih m hmlt hmpos hexpM hcardM
      refine ⟨hsplit ▸ ((zmodProdCrt (p ^ k) m hco).symm.trans ((eP.prodCongr eM).trans eCrt))⟩

/-! ## ★出典の紐付け(`.src`) -/

def addEquiv_zmod_sq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(有限アーベル群の構造定理——捩れが (Z/n)²)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
