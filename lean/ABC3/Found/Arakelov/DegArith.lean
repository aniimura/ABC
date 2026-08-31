/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.ArchDegSmul
import ABC3.Found.Arakelov.APicToSheaf
import ABC3.Found.Arakelov.PicAssoc
import ABC3.Meta.Claim

/-!
# **`deg_F` の古典的な定義**と、その**切断の取り方への依らなさ**（段 D、`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

## ★★★★★★★★★★これは何か

    `deg_F(L̄) = ( log #(Gamma(L)/O_F·s) − Sum_sigma log |s|_sigma ) / [F:Q]`

本ファイルはこの式を `degArith` として**書き**、
`s` を `c·s`（`c ∈ O_F`、`c ≠ 0`）に取り替えても**変わらない**ことを示す。

| 側 | ずれ | 場所 |
|---|---|---|
| 有限側 | `+ log |N(c)|` | `degFin_smul_norm`（§9-774） |
| アルキメデス側 | `− log |N(c)| / [F:Q]` | `archDeg_smul`（§9-775） |
| ★**和** | `0` | ★**本ファイル**（§9-776） |

★★★これが**積公式**の内容である——原文が `deg_F` を「切断の取り方に依らない」量として
使えるのは、この相殺が起きるからである。

## ★★★★★二つの世界をつなぐ蝶番

`degFin` は**層化した加群** `Gamma(sheafify L)` の上に、
`archDeg` は**前層の切断** `L.sheaf.obj (op ⊤)` の上にある。
★`gammaSheafifyM`（層化の単位を `⊤` で評価したもの）がその橋である。
★★`gammaSheafifyM_smul` が「係数は `Gamma-Spec` 同型で移る」ことを言う
——これで両者の `c` が同じ `c` になる。

## ★残っている段（明示）

★★本ファイルが与えるのは「**大域函数倍で不変**」である。
★★★**任意の二つの非零切断が大域函数倍で結ばれること**（前層 `Gamma` の階数 1 性）は
まだ入っていない——層化した側では `exists_smul_mem_span_gamma`（段 B）が
それを言うが、**前層の切断の水準では層化の単位の単射性が要る**。
局所自明性だけからは前層の分離性は従わないので、これは独立の段である。
★★★★`degArith_congr_of_smul_eq` は「共通の倍元をもつ二つの切断」については
すでに結論を与えている。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite NumberField

/-! ## ★★★★★★層化の単位を `Γ` で評価する -/

/-- ★★★★★★**前層の大域切断を層化した加群へ送る**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★`degFin` は層化した側（`moduleSpecΓFunctor`）に、
`archDeg` は前層の側（`L.sheaf.obj (op ⊤)`）にある。★★これがその橋である。 -/
noncomputable def gammaSheafifyM (R : CommRingCat.{0}) (L : AInv (Spec R))
    (s : (L.carrier.sheaf.obj (op ⊤) : Type)) :
    (AlgebraicGeometry.moduleSpecΓFunctor.obj (AInv.toInvSheaf L).carrier : Type) :=
  ((sheafifyUnit (Spec R) L.carrier.sheaf).app (op ⊤)).hom s

/-- ★★★★★★★**橋は係数と両立する**——前層の側の `Γ(X,⊤)` 倍は、
層化した側では `Γ-Spec` 同型で移した `R` 倍になる。

★機構は `PresheafOfModules` の射が `Γ(X,⊤)`-線形であることと、
`R`-作用が `ΓSpecIso.inv` に沿った係数制限であること（`rfl`）である。 -/
theorem gammaSheafifyM_smul (R : CommRingCat.{0}) (L : AInv (Spec R))
    (c : (Γ(Spec R, (⊤ : (Spec R).Opens)) : Type))
    (s : (L.carrier.sheaf.obj (op ⊤) : Type)) :
    gammaSheafifyM R L (c • s) = ((Scheme.ΓSpecIso R).hom.hom c) • gammaSheafifyM R L s := by
  have hc : (Scheme.ΓSpecIso R).inv.hom ((Scheme.ΓSpecIso R).hom.hom c) = c :=
    congrArg (fun (m : _ ⟶ _) => CommRingCat.Hom.hom m c) (Scheme.ΓSpecIso R).hom_inv_id
  have h : ((sheafifyUnit (Spec R) L.carrier.sheaf).app (op ⊤)).hom (c • s)
      = c • ((sheafifyUnit (Spec R) L.carrier.sheaf).app (op ⊤)).hom s :=
    ((sheafifyUnit (Spec R) L.carrier.sheaf).app (op ⊤)).hom.map_smul c s
  show ((sheafifyUnit (Spec R) L.carrier.sheaf).app (op ⊤)).hom (c • s) = _
  rw [h]
  exact congrArg (fun r => r • (((sheafifyUnit (Spec R) L.carrier.sheaf).app (op ⊤)).hom s)) hc.symm

/-! ## ★★★★★★★★★★`deg_F` の古典的な定義 -/

/-- ★**`𝓞_F` の非零元による商は有限である**（`degFin_smul_norm` の欄を埋める）。 -/
theorem finite_quotient_span (F : Type) [Field F] [NumberField F]
    (r : ((CommRingCat.of (𝓞 F)) : Type)) (hr : r ≠ 0) :
    Finite (((CommRingCat.of (𝓞 F)) : Type) ⧸ (Ideal.span {r} : Ideal (𝓞 F))) := by
  have h : (Ideal.span {r} : Ideal (𝓞 F)) ≠ ⊥ := by
    simpa [Ideal.span_singleton_eq_bot] using hr
  exact Ideal.finiteQuotientOfFreeOfNeBot _ h

/-- ★★★★★★★★★★**`deg_F(L̄)`（切断 `s` を選んだ形）**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

    `degArith L̄ s = ( log #(Γ(L)/𝓞_F·s) − Σ_σ log |s|_σ ) / [F:ℚ]`

★★第 1 項が `degFin`（段 B）、第 2 項が `archDeg`（段 C）である。
★★★正規化 `/[F:ℚ]` は原文の `deg_F` の規約（底変換で不変になる形）に合わせてある。 -/
noncomputable def degArith (F : Type) [Field F] [NumberField F]
    (L : AInv (Spec (CommRingCat.of (𝓞 F))))
    (s : (L.carrier.sheaf.obj (op ⊤) : Type)) : ℝ :=
  degFin (CommRingCat.of (𝓞 F)) L (gammaSheafifyM (CommRingCat.of (𝓞 F)) L s)
      / (Module.finrank ℚ F : ℝ)
    + archDeg F L.carrier s

/-! ## ★★★★★★★★★★段 D —— 切断を大域函数倍しても変わらない -/

/-- ★★★★★★★★★★**段 D**——`deg_F` は切断の大域函数倍で不変である。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★有限側は `+ log |N(c)|`（`degFin_smul_norm`、§9-774）、
アルキメデス側は `− log |N(c)|/[F:ℚ]`（`archDeg_smul`、§9-775）動く。
★★正規化を揃えると**ちょうど打ち消す**——これが**積公式**である。

★★★これで台帳 `arakelov-degF-finite-places` の段 D の中核が入った。 -/
theorem degArith_smul (F : Type) [Field F] [NumberField F]
    (L : AInv (Spec (CommRingCat.of (𝓞 F))))
    (c : (Γ(Spec (CommRingCat.of (𝓞 F)),
      (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type))
    (s : (L.carrier.sheaf.obj (op ⊤) : Type))
    (hc : (Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).hom.hom c ≠ 0)
    (hgs : gammaSheafifyM (CommRingCat.of (𝓞 F)) L s ≠ 0)
    (hs : ∀ σ : F →+* ℂ, L.carrier.norm s (embSpecPoint F σ) ≠ 0) :
    degArith F L (c • s) = degArith F L s := by
  have hfin := degFin_smul_norm (CommRingCat.of (𝓞 F)) L
    (gammaSheafifyM (CommRingCat.of (𝓞 F)) L s) hgs
    ((Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).hom.hom c) hc (finite_quotient_span F)
  have harch := archDeg_smul F L.carrier c s hc hs
  show degFin (CommRingCat.of (𝓞 F)) L
      (gammaSheafifyM (CommRingCat.of (𝓞 F)) L (c • s)) / _ + archDeg F L.carrier (c • s) = _
  rw [gammaSheafifyM_smul, hfin, harch]
  show _ = degFin (CommRingCat.of (𝓞 F)) L
      (gammaSheafifyM (CommRingCat.of (𝓞 F)) L s) / _ + archDeg F L.carrier s
  ring

/-- ★★★★★★★★**共通の倍元をもつ二つの切断は同じ次数を与える**。

★これが「`deg_F` は切断の取り方に依らない」の**使える形**である
——`a·s = b·t` なら `deg(s) = deg(t)`。
★★層化した側では `exists_smul_mem_song_gamma`（段 B）が
任意の二つの非零切断についてそのような `a`, `b` を与える。 -/
theorem degArith_congr_of_smul_eq (F : Type) [Field F] [NumberField F]
    (L : AInv (Spec (CommRingCat.of (𝓞 F))))
    (a b : (Γ(Spec (CommRingCat.of (𝓞 F)),
      (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type))
    (s t : (L.carrier.sheaf.obj (op ⊤) : Type))
    (ha : (Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).hom.hom a ≠ 0)
    (hb : (Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).hom.hom b ≠ 0)
    (hgs : gammaSheafifyM (CommRingCat.of (𝓞 F)) L s ≠ 0)
    (hgt : gammaSheafifyM (CommRingCat.of (𝓞 F)) L t ≠ 0)
    (hs : ∀ σ : F →+* ℂ, L.carrier.norm s (embSpecPoint F σ) ≠ 0)
    (ht : ∀ σ : F →+* ℂ, L.carrier.norm t (embSpecPoint F σ) ≠ 0)
    (heq : a • s = b • t) :
    degArith F L s = degArith F L t := by
  have h1 := degArith_smul F L a s ha hgs hs
  have h2 := degArith_smul F L b t hb hgt ht
  rw [← h1, ← h2, heq]

/-! ### ★出典の紐付け(`.src`) -/

def gammaSheafifyM.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(前層の大域切断を層化した加群へ送る橋)",
    sectionId := "genell-def-1-1-ii" }

def degArith.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(deg_F の古典的な定義——有限部分＋アルキメデス部分)",
    sectionId := "genell-def-1-1-ii" }

def degArith_smul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(段 D——切断の大域函数倍で deg_F は不変＝積公式)",
    sectionId := "genell-def-1-1-ii" }

def degArith_smul.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "degFin_smul_norm(有限側のずれ = + log |N(c)|、§9-774)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.degFin_smul_norm") 4,
    .citation "[ABC3]" "archDeg_smul(アルキメデス側のずれ = − log |N(c)|/[F:Q]、§9-775)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.archDeg_smul") 4,
    .citation "[mathlib]" "Ideal.finiteQuotientOfFreeOfNeBot(非零イデアルによる商は有限)"
      (.inMathlib "Ideal.finiteQuotientOfFreeOfNeBot") 4,
    .implicitStep
      ("★★本定理が与えるのは「大域函数倍で不変」である。" ++
       "★★★任意の二つの非零切断が大域函数倍で結ばれること(前層 Γ の階数 1 性)は" ++
       "まだ入っていない——層化した側では exists_smul_mem_span_gamma(段 B)が" ++
       "それを言うが、前層の切断の水準では層化の単位の単射性が要る。" ++
       "局所自明性だけからは前層の分離性は従わないので、これは独立の段である") 4 ]

end ABC3.Found.Arakelov
