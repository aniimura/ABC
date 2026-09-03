/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateSpecialize
import ABC3.Found.GenEll.Velu
import ABC3.Found.GenEll.JScale
import ABC3.Found.GaloisRep.TateEquation
import ABC3.Meta.Claim
import ABC3.Found.GaloisRep.TateVelu.Definition33

/-!
# TateVelu —— `[GenEll] Lemma 3.5` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve ABC3.Found.GenEll
variable {R : Type} [CommRing R] {I : Ideal R}

def veluGx_tateCurveAt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Tate 曲線では g^x_Q = 3x² + a₄ − y。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluV2_tateCurveAt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(v_Q = g^x_Q の Tate 形。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluGy_tateCurveAt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Tate 曲線では g^y_Q = −2y − x。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluU_tateCurveAt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(u_Q = (2y + x)² の Tate 形。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
